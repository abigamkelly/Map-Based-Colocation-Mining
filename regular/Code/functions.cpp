#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <chrono>
#include <cmath>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/range/algorithm/sort.hpp>

extern "C" {
    namespace fs = std::filesystem;
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;
    typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;  // cartesian point
    typedef std::pair<point_t, size_t> rtree_value_t;  // rtree value type of point w index
    typedef bg::model::box<point_t> Box;
    typedef bg::model::polygon<point_t> Polygon;
    
    struct Feature {
        int feature_type;
        int id;
        int pos;
        double latitude;
        double longitude;

        bool operator<(const Feature& other) const {
            return feature_type < other.feature_type;
        }
    };
    
    struct FeatureInformation {
        int count;
        int start;
        int end;
    };
    
    class Region {
        public:
        std::map<int, FeatureInformation> featureInfo;
        std::map<int, std::vector<int>> starNeighbors;
        std::vector<std::vector<int>> size2_patterns;
        std::map<std::vector<int>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<int>, std::map<int, std::set<int>>> hashmap;
    };
    
    Region region;
    float distance_threshold = 10;
    float prevalence_threshold = 0.5;
    
    float calc_distance(const Feature& f1, const Feature& f2) {
        return std::sqrt(std::pow(f1.latitude - f2.latitude, 2) + 
                         std::pow(f1.longitude - f2.longitude, 2));
    }
    
    std::vector<Feature> read_data(const std::string& folderPath) {
        std::vector<Feature> features;
        int region_count = 0;
        int row_count = 0;
        int total_count = 0;


        std::ifstream file(folderPath);
        if (!file.is_open()) {
            std::cerr << "Error opening file!" << std::endl;
            exit(0);
        }

        std::string line;
        std::getline(file, line);

        while (std::getline(file, line)) {
            Feature tempFeature;
            tempFeature.pos = row_count;
            row_count += 1;
            total_count += 1;
            std::stringstream ss(line);
            std::string value;
            int colIndex = 0;


            while (std::getline(ss, value, ',')) {
                if (colIndex == 0) {
                    tempFeature.feature_type = std::stoi(value); 
                }
                else if (colIndex == 1) {
                    tempFeature.latitude = std::stod(value);
                }
                else if (colIndex == 2) {
                    tempFeature.longitude = std::stod(value);
                }
                colIndex++;
            }
                        
            features.push_back(tempFeature);
        }
        file.close();

        return features;
    }
    
    std::map<int, FeatureInformation> calc_featureInfo(std::vector<Feature> features) {
        std::map<int, FeatureInformation> featureInfo;
        int count = 0;
        int start_row_id = 0;
        int prev_feature = -1;

        for(int i = 0; i < features.size(); i++) {
            int feature = features[i].feature_type;
            if(feature != prev_feature) {
                if(prev_feature != -1) {
                    featureInfo[prev_feature] = {count, start_row_id, i - 1};
                    //std::cout << prev_feature << " " << count << " " << start_row_id << " " << i - 1 << std::endl;
                }
                count = 1;
                start_row_id = i;
                prev_feature = feature;
            } else {
                count++; 
            }
        }
        if(prev_feature != -1) {
            featureInfo[prev_feature] = {count + 1, start_row_id, start_row_id + count - 1};
            //std::cout << prev_feature << " " << count + 1 << " " << start_row_id << " " << start_row_id + count - 1 << std::endl;
        }
        if(prev_feature != -1 && count == 1) {
            featureInfo[prev_feature] = {count, start_row_id, start_row_id};
            //std::cout << prev_feature << " " << count << " " << start_row_id << " " << start_row_id << std::endl;
        }  
        return featureInfo;
    }
    
    
    std::map<int, std::vector<int>> calc_starNeighbors(std::vector<Feature> features) {
        std::map<int, std::vector<int>> starNeighbors;

        std::set<int> feature_types;
        for (const auto& feature : features) {
            feature_types.insert(feature.feature_type);
        }

        for(int current_feature : feature_types) {
            std::vector<Feature> current_feature_coords;
            std::vector<Feature> non_current_coords;

            for (const auto& feature : features) {
                // get coordinates of current feature type
                if (feature.feature_type == current_feature) {
                    current_feature_coords.push_back(feature);
                } else {  // get coordinates of other feature types
                    non_current_coords.push_back(feature);
                }
            }

            // create rtree
            bgi::rtree<rtree_value_t, bgi::quadratic<16>> rtree;

            // insert points of non-current feature type into rtree
            for (size_t i = 0; i < non_current_coords.size(); ++i) {
                point_t pt(non_current_coords[i].latitude, non_current_coords[i].longitude);
                rtree.insert(std::make_pair(pt, non_current_coords[i].pos));  // Insert point with index i
            }

            for(size_t i = 0; i < current_feature_coords.size(); ++i) {
                // Define the center point
                point_t center(current_feature_coords[i].latitude, current_feature_coords[i].longitude);

                Box query_box(point_t(current_feature_coords[i].latitude - distance_threshold, 
                                    current_feature_coords[i].longitude - distance_threshold),
                              point_t(current_feature_coords[i].latitude + distance_threshold, 
                                    current_feature_coords[i].longitude + distance_threshold));
                std::vector<rtree_value_t> nearby_points;
                rtree.query(bgi::intersects(query_box), std::back_inserter(nearby_points));

                std::vector<int> points_to_add;         
                for (const auto& neighbor : nearby_points) {
                    int neighbor_id = neighbor.second;
                    const Feature& neighbor_feature = *std::find_if(non_current_coords.begin(), non_current_coords.end(),
                        [&neighbor_id](const Feature& f) { return f.pos == neighbor_id; });

                    // Ensure the neighbor is farther, of a different feature type, and within the distance threshold
                    if (neighbor_id > current_feature_coords[i].pos && calc_distance(current_feature_coords[i], neighbor_feature) <= distance_threshold) {
                        points_to_add.push_back(neighbor_id);
                    }
                }

                std::sort(points_to_add.begin(), points_to_add.end());
                starNeighbors[current_feature_coords[i].pos] = points_to_add; 
            }
        }
        return starNeighbors;
    }
    
    std::vector<std::vector<int>> generate_size2_combos(std::map<int, FeatureInformation> featureInfo) {
        std::vector<int> features;
        for (const auto &entry : featureInfo) {
            features.push_back(entry.first);
        }

        std::vector<std::vector<int>> size2_candidatePatterns;
        for (size_t i = 0; i < features.size(); ++i) {
            for (size_t j = i + 1; j < features.size(); ++j) {
                size2_candidatePatterns.push_back({features[i], features[j]});
            }
        }

        return size2_candidatePatterns;
    }
    
    std::vector<int> findNeighborsInRange(const std::vector<int>& arr, int x, int y) {
        auto start_it = std::lower_bound(arr.begin(), arr.end(), x);
        auto end_it = std::upper_bound(arr.begin(), arr.end(), y);
        return std::vector<int>(start_it, end_it);
    }
    
    
    std::vector<std::vector<int>> degree2Processing(std::vector<std::vector<int>> candidatePatterns) {   
        std::vector<std::vector<int>> size2_patterns;

        for (const auto& coloc : candidatePatterns) {
            int first_feature = coloc[0];
            int second_feature = coloc[1];

            // determine first feature start and end
            auto feature1 = region.featureInfo.find(first_feature);
            FeatureInformation values1 = feature1->second;
            int first_feature_start = values1.start;
            int first_feature_end = values1.end;
            // determine second feature start and end
            auto feature2 = region.featureInfo.find(second_feature);
            FeatureInformation values2 = feature2->second;
            int second_feature_start = values2.start;
            int second_feature_end = values2.end;

            std::vector<int> coloc_key = {first_feature, second_feature};
            region.instance_table[coloc_key] = {};
            region.hashmap[coloc_key] = {};
            region.hashmap[coloc_key][first_feature] = {};
            region.hashmap[coloc_key][second_feature] = {};

            for (int index = first_feature_start; index <= first_feature_end; index++) {
                auto star_neighbor_it = region.starNeighbors.find(index);
                std::vector<int> neighbors = findNeighborsInRange(star_neighbor_it->second, 
                                                                  second_feature_start,
                                                                  second_feature_end);
                if (!neighbors.empty()) {
                    std::vector<int> index_tuple = {index};
                    region.instance_table[coloc_key][index_tuple] = neighbors;
                    region.hashmap[coloc_key][first_feature].insert(index);
                    for (int neighbor : neighbors) {
                        region.hashmap[coloc_key][second_feature].insert(neighbor);
                    }
                }
            }

            double pr_first_feature = static_cast<double>(region.hashmap[coloc_key][first_feature].size()) /
                region.featureInfo[first_feature].count;
            double pr_second_feature = static_cast<double>(region.hashmap[coloc_key][second_feature].size()) /
                region.featureInfo[second_feature].count;

            double PI = 0.0;

            if (pr_first_feature < pr_second_feature) {
                PI = pr_first_feature;
            } else {
                PI = pr_second_feature;
            }

            if (PI >= prevalence_threshold) {
                //std::cout << hashmap[coloc_key][first_feature].size() << " / " << featureInfo[first_feature].count << " , " << hashmap[coloc_key][second_feature].size() << " / " << featureInfo[second_feature].count << std::endl;
                   size2_patterns.push_back(coloc_key);
            }
            else {
                region.instance_table.erase(coloc_key);
            }
        }

        std::cout << "Degree 2 Prevalent Patterns:" << std::endl;
        for (auto i : size2_patterns) {
            std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
        }
        return size2_patterns;
    }
    
    
    void generateCombinations(const std::vector<int>& features, int degree, 
                              std::vector<std::vector<int>>& result, 
                              std::vector<int>& current, int start) {
        if (current.size() == degree) {
            result.push_back(current);
            return;
        }
        for (int i = start; i < features.size(); ++i) {
            current.push_back(features[i]);
            generateCombinations(features, degree, result, current, i + 1);
            current.pop_back();
        }
    }

    // Function to check if all (degree-1)-subpatterns of a pattern are in the prevalent patterns
    bool allSubpatternsInPrevalent(const std::vector<int>& pattern, 
                                   const std::set<std::vector<int>>& prevalentPatterns, 
                                   int degree) {
        std::vector<std::vector<int>> subpatterns;
        std::vector<int> current;
        generateCombinations(pattern, degree - 1, subpatterns, current, 0);

        for (const auto& subpattern : subpatterns) {
            if (prevalentPatterns.find(subpattern) == prevalentPatterns.end()) {
                return false;
            }
        }
        return true;
    }

    std::vector<std::vector<int>> getCandidatePatterns(const std::vector<std::vector<int>>& prevalentPattern, 
                                                       int degree) {
        std::set<std::vector<int>> prevalentPatterns(prevalentPattern.begin(), prevalentPattern.end());
        // Extract features from the keys of featureInfo
        std::vector<int> features;
        for (const auto& pair : region.featureInfo) {
            features.push_back(pair.first);
        }

        std::vector<std::vector<int>> _candidatePatterns;
        std::vector<std::vector<int>> _patterns;
        std::vector<int> current;
        generateCombinations(features, degree, _patterns, current, 0);

        for (const auto& pattern : _patterns) {
            if (allSubpatternsInPrevalent(pattern, prevalentPatterns, degree)) {
                _candidatePatterns.push_back(pattern);
            }
        }
        return _candidatePatterns;
    }
    
    
    std::vector<std::vector<int>> colocationGeneral(std::vector<std::vector<int>> candidatePatterns, int degree) {
        unsigned long new_key_count = 0;
        unsigned long inner_list_count = 0;
        std::vector<std::vector<int>> prevalent;

        for (const auto& currentPattern : candidatePatterns) {
            std::vector<int> basePattern;
            for (int j = 0; j < degree - 1; j++) {
                basePattern.push_back(currentPattern[j]);
            }
            int lastFeature = currentPattern[degree - 1];
            // Add entries to instance_table and hashmap
            region.instance_table[currentPattern] = {};
            region.hashmap[currentPattern] = {};

            // Initialize sets for each element in currentPattern in hashmap
            for (const auto& f : currentPattern) {
                region.hashmap[currentPattern][f] = {};
            }

            auto colocTableIt = region.instance_table.find(basePattern);
            std::map<std::vector<int>, std::vector<int>>& colocTable = colocTableIt->second;
            
            
            for (const auto& entry : colocTable) {
                const std::vector<int>& key = entry.first;
                std::set<int> commonLastNeighbors;
                for (int instanceID : key) {
                    auto star_neighbor_it = region.starNeighbors.find(instanceID);
                    if (commonLastNeighbors.empty()) {
                        std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second,
                                                                     region.featureInfo[lastFeature].start,
                                                                     region.featureInfo[lastFeature].end);
                        commonLastNeighbors.insert(temp_vector.begin(), temp_vector.end());
                    } else {
                        std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second, 
                                                                              region.featureInfo[lastFeature].start,
                                                                              region.featureInfo[lastFeature].end);
                        std::set<int> temp_commonLastNeighbors(temp_vector.begin(), temp_vector.end());
                        std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                              temp_commonLastNeighbors.begin(), temp_commonLastNeighbors.end(),
                                              std::inserter(commonLastNeighbors, commonLastNeighbors.begin()));
                    }
                }
                for (int n : colocTable[key]) {
                    auto star_neighbor_it = region.starNeighbors.find(n);
                    std::vector<int> temp_vect = findNeighborsInRange(star_neighbor_it->second,
                                                                      region.featureInfo[lastFeature].start,
                                                                      region.featureInfo[lastFeature].end);
                    std::set<int> temp_neighbors(temp_vect.begin(), temp_vect.end());
                    std::set<int> neighbors;
                    std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                          temp_neighbors.begin(), temp_neighbors.end(),
                                          std::inserter(neighbors, neighbors.begin()));

                    if (!neighbors.empty()) {
                        std::vector<int> new_key = key;
                        new_key.push_back(n);
                        
                        
                        if (degree < 4) {
                            std::vector<int> intersectionVec(neighbors.begin(), neighbors.end());
                            region.instance_table[currentPattern][new_key] = intersectionVec;
                        }
                        else {
                            new_key_count += 1;
                            inner_list_count += neighbors.size();
                        }

                        for (size_t k = 0; k < new_key.size(); ++k) {
                            region.hashmap[currentPattern][currentPattern[k]].insert(new_key[k]);
                        }
                        region.hashmap[currentPattern][lastFeature].insert(neighbors.begin(), neighbors.end());
                    }
                }
            }

            std::vector<double> pr;
            for (int m = 0; m < degree; ++m) {
                int f = currentPattern[m];
                double ratio = static_cast<double>(region.hashmap[currentPattern][f].size()) 
                    / region.featureInfo[f].count;
                pr.push_back(ratio);
            }
            double PI = *std::min_element(pr.begin(), pr.end());
            if (PI >= prevalence_threshold) {
                prevalent.push_back(currentPattern);
            }
            else {
                region.instance_table.erase(currentPattern);
            }
        }
        
        if (degree >= 4) {
            std::cout << "NEW_KEY_COUNT: " << new_key_count << " INNER_LIST_COUNT: " << inner_list_count << std::endl;
        }

        if(!prevalent.empty()) {
            std::cout << "Degree " << degree << " Prevalent Patterns: " << std::endl;

            for (auto& patternVec : prevalent) {
                std::cout << "(";
                for (size_t i = 0; i < patternVec.size(); i++) {
                    std::cout << patternVec[i];
                    if (i < patternVec.size() - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << ")" << std::endl;
            }
        }
        return prevalent;
    }
    
    
    int main() {
        std::string input_file = "/home/amk7r/colocation_mining/TestCase2/TestCase2_1.csv";
        std::vector<Feature> features = read_data(input_file);
        
        auto start = std::chrono::high_resolution_clock::now();
        
        region.featureInfo = calc_featureInfo(features);
        std::cout << "featureInfo\n";
        region.starNeighbors = calc_starNeighbors(features);
        std::cout << "starNeighbors\n";
        
        std::vector<std::vector<int>> size2_candidatePatterns = generate_size2_combos(region.featureInfo);
        region.size2_patterns = degree2Processing(size2_candidatePatterns);
        
        int key_count = 0;
        int inner_key_count = 0;
        int inner_value_count = 0;

        // Iterate through the instance_table
        for (const auto& outer_pair : region.instance_table) {
            const std::vector<int>& key = outer_pair.first;
            const std::map<std::vector<int>, std::vector<int>>& inner_map = outer_pair.second;

            // Check if the size of the key vector is 2
            if (key.size() == 2) {
                key_count++;

                // Iterate through the inner map
                for (const auto& inner_pair : inner_map) {
                    const std::vector<int>& inner_key = inner_pair.first;
                    const std::vector<int>& inner_value = inner_pair.second;

                    inner_key_count++;
                    inner_value_count += inner_value.size();
                }
            }
        }
        
        std::cout << "Degree: " << 2 << ", Number of keys: " << key_count << ", Number of inner keys: " << inner_key_count << ", Number of inner values: " << inner_value_count << std::endl;
        
        
        
        
        int degree = 3;
        std::vector<std::vector<int>> candidatePatterns = getCandidatePatterns(region.size2_patterns, degree);
        
        while (!candidatePatterns.empty()) {
            std::vector<std::vector<int>> prevalent_patterns = colocationGeneral(candidatePatterns, degree);
            
            unsigned long key_count = 0;
            unsigned long inner_key_count = 0;
            unsigned long inner_value_count = 0;

            // Iterate through the instance_table
            for (const auto& outer_pair : region.instance_table) {
                const std::vector<int>& key = outer_pair.first;
                const std::map<std::vector<int>, std::vector<int>>& inner_map = outer_pair.second;

                // Check if the size of the key vector is 2
                if (key.size() == degree) {
                    key_count++;

                    // Iterate through the inner map
                    for (const auto& inner_pair : inner_map) {
                        const std::vector<int>& inner_key = inner_pair.first;
                        const std::vector<int>& inner_value = inner_pair.second;

                        inner_key_count++;
                        inner_value_count += inner_value.size();
                    }
                }
            }
            
            std::cout << "Degree: " << degree << ", Number of keys: " << key_count << ", Number of inner keys: " << inner_key_count << ", Number of inner values: " << inner_value_count << std::endl;
            
            
            
            degree += 1;
            if (prevalent_patterns.size() == 0) {
                break;
            }

            if(degree == 5) {
                break;
            }
            candidatePatterns = getCandidatePatterns(prevalent_patterns, degree);
        }
        
        
        auto end = std::chrono::high_resolution_clock::now();

        // Calculate duration in different units
        auto duration_sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        auto duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

        std::cout << "Time taken: " << duration_sec.count() << " seconds\n";
        std::cout << "Time taken: " << duration_ms.count() << " milliseconds\n";
        std::cout << "Time taken: " << duration_ns.count() << " nanoseconds\n";
        
        
        return 0;
    }
    
    
}
