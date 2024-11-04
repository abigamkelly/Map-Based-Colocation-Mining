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
    
    class Subregion {
        public:
        int subregion_id;
        std::vector<Feature> data;
        std::map<int, FeatureInformation> featureInfo;
        std::map<int, std::vector<int>> starNeighbors;
        std::vector<std::vector<int>> size2_patterns;
        std::map<std::vector<int>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<int>, std::map<int, std::set<int>>> hashmap;
        
        Subregion(int number) { 
            subregion_id = number; 
        }
    };

    class Border {
        public:
        int border_id;
        std::map<int, FeatureInformation> featureInfo;
        std::map<int, std::vector<int>> starNeighbors;
        std::vector<std::vector<int>> size2_patterns;
        std::map<std::vector<int>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<int>, std::map<int, std::set<int>>> hashmap;
        
        Border(int number) {
            border_id = number;
        }
    };
    
    class Region {
        public:
        std::vector<std::vector<int>> size2_patterns;
        std::map<std::vector<int>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<int>, std::map<int, std::set<int>>> hashmap;
    };
    
    float distance_threshold = 24.53;
    float prevalence_threshold = 0.55;
    std::vector<Subregion> subregions;
    std::vector<Border> borders;
    Region region;
    
    int main() {
        return 0;
    }
    
    float read_distance(const std::string& distance_file) {
        std::ifstream infile(distance_file);
        if (!infile.is_open()) {
            std::cerr << "Error: Could not open file: " << distance_file << std::endl;
            return 1;
        }
        infile >> distance_threshold;
        infile.close();
        return distance_threshold;
    }


    std::vector<Feature> read_data(const std::string& folderPath, std::vector<int> &offsets,
                                  std::vector<int> &indices) {
        std::vector<Feature> features;
        int region_count = 0;
        int row_count = 0;
        int total_count = 0;

        for (const auto& entry : fs::directory_iterator(folderPath)) {
            if (entry.path().extension() == ".csv") {
                region_count += 1;
                row_count = 0;
                indices.push_back(row_count);
                std::ifstream file(entry.path());
                if (!file.is_open()) {
                    std::cerr << "Error opening file: " << entry.path() << std::endl;
                    continue;
                }

                std::string line;
                std::getline(file, line);

                while (std::getline(file, line)) {
                    Feature tempFeature;
                    tempFeature.id = total_count;
                    row_count += 1;
                    total_count += 1;
                    std::stringstream ss(line);
                    std::string value;
                    int colIndex = 0;


                    while (std::getline(ss, value, ',')) {
                        if (colIndex == 0) {
                            tempFeature.feature_type = std::stoi(value); 
                        }
                        else if (colIndex == 3) {
                            tempFeature.latitude = std::stod(value);
                        }
                        else if (colIndex == 4) {
                            tempFeature.longitude = std::stod(value);
                        }
                        colIndex++;
                    }
                    features.push_back(tempFeature);
                }
                indices.push_back(row_count - 1);
                file.close();
                offsets.push_back(row_count);
            }
        }
        return features;
    }
    
    std::map<int, FeatureInformation> calc_featureInfo(std::vector<Feature> features,
                                                   int start_idx, int end_idx, int offset) {
        std::map<int, FeatureInformation> featureInfo;
        int count = 0;
        int start_row_id = 0;
        int prev_feature = -1;

        for(int i = start_idx + offset; i < end_idx + offset; i++) {
            int feature = features[i].feature_type;
            if(feature != prev_feature) {
                if(prev_feature != -1) {
                    featureInfo[prev_feature] = {count, start_row_id, i - 1};
                }
                count = 1;
                start_row_id = i;
                prev_feature = feature;
            } else {
                count++; 
            }
        }
        if(prev_feature != -1) {
            featureInfo[prev_feature] = {count + 1, start_row_id, start_idx + end_idx + offset};
        }
        if(prev_feature != -1 && count == 1) {
            featureInfo[prev_feature] = {count, start_row_id, start_row_id};
        }  
        return featureInfo;
    }
    
    float calc_distance(const Feature& f1, const Feature& f2) {
        return std::sqrt(std::pow(f1.latitude - f2.latitude, 2) + 
                         std::pow(f1.longitude - f2.longitude, 2));
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
                rtree.insert(std::make_pair(pt, non_current_coords[i].id));  // Insert point with index i
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
                        [&neighbor_id](const Feature& f) { return f.id == neighbor_id; });

                    // Ensure the neighbor is farther, of a different feature type, and within the distance threshold
                    if (neighbor_id > current_feature_coords[i].id && calc_distance(current_feature_coords[i], neighbor_feature) <= distance_threshold) {
                        points_to_add.push_back(neighbor_id);
                    }
                }

                std::sort(points_to_add.begin(), points_to_add.end());
                starNeighbors[current_feature_coords[i].id] = points_to_add; 
            }
        }
        return starNeighbors;
    }
    
    
    std::vector<int> findNeighborsInRange(const std::vector<int>& arr, int x, int y) {
        auto start_it = std::lower_bound(arr.begin(), arr.end(), x);
        auto end_it = std::upper_bound(arr.begin(), arr.end(), y);
        return std::vector<int>(start_it, end_it);
    }
    
    
    std::vector<std::vector<int>> colocationGeneral(
        std::vector<std::vector<int>> candidatePatterns, double prevalence_threshold, 
        int degree, int number_subregions, 
        std::map<std::vector<int>, std::map<std::vector<int>, std::vector<int>>>& instance_table, 
        std::map<std::vector<int>, std::map<int, std::set<int>>>& hashmap, 
        std::map<int, std::vector<int>> starNeighbors, std::map<int, FeatureInformation> featureInfo) {
        
        std::vector<std::vector<int>> prevalent;

        for (const auto& currentPattern : candidatePatterns) {
            std::vector<int> basePattern;
            for (int j = 0; j < degree - 1; j++) {
                basePattern.push_back(currentPattern[j]);
            }
            int lastFeature = currentPattern[degree - 1];
            // Add entries to instance_table and hashmap
            instance_table[currentPattern] = {};
            hashmap[currentPattern] = {};

            // Initialize sets for each element in currentPattern in hashmap
            for (const auto& f : currentPattern) {
                hashmap[currentPattern][f] = {};
            }

            auto colocTableIt = instance_table.find(basePattern);
            std::map<std::vector<int>, std::vector<int>>& colocTable = colocTableIt->second;
            
            
            for (const auto& entry : colocTable) {
                const std::vector<int>& key = entry.first;
                std::set<int> commonLastNeighbors;
                for (int instanceID : key) {
                    auto star_neighbor_it = starNeighbors.find(instanceID);
                    if (commonLastNeighbors.empty()) {
                        std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second,
                                                                     featureInfo[lastFeature].start,
                                                                     featureInfo[lastFeature].end);
                        commonLastNeighbors.insert(temp_vector.begin(), temp_vector.end());
                    } else {
                        std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second, 
                                                                              featureInfo[lastFeature].start,
                                                                              featureInfo[lastFeature].end);
                        std::set<int> temp_commonLastNeighbors(temp_vector.begin(), temp_vector.end());
                        std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                              temp_commonLastNeighbors.begin(), temp_commonLastNeighbors.end(),
                                              std::inserter(commonLastNeighbors, commonLastNeighbors.begin()));
                    }
                }
                for (int n : colocTable[key]) {
                    auto star_neighbor_it = starNeighbors.find(n);
                    std::vector<int> temp_vect = findNeighborsInRange(star_neighbor_it->second,
                                                                      featureInfo[lastFeature].start,
                                                                      featureInfo[lastFeature].end);
                    std::set<int> temp_neighbors(temp_vect.begin(), temp_vect.end());
                    std::set<int> neighbors;
                    std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                          temp_neighbors.begin(), temp_neighbors.end(),
                                          std::inserter(neighbors, neighbors.begin()));

                    if (!neighbors.empty()) {
                        std::vector<int> new_key = key;
                        new_key.push_back(n);
                        std::vector<int> intersectionVec(neighbors.begin(), neighbors.end());
                        instance_table[currentPattern][new_key] = intersectionVec;

                        for (size_t k = 0; k < new_key.size(); ++k) {
                            hashmap[currentPattern][currentPattern[k]].insert(new_key[k]);
                        }
                        hashmap[currentPattern][lastFeature].insert(neighbors.begin(), neighbors.end());
                    }
                }
            }

            std::vector<double> pr;
            for (int m = 0; m < degree; ++m) {
                int f = currentPattern[m];
                double ratio = static_cast<double>(hashmap[currentPattern][f].size()) 
                    / featureInfo[f].count;
                pr.push_back(ratio);
            }
            double PI = *std::min_element(pr.begin(), pr.end());
            if (PI >= prevalence_threshold) {
                prevalent.push_back(currentPattern);
            }
        }

        if(!prevalent.empty()) {
            std::cout << "Degree " << degree << " Prevalent Patterns for Sub-Region " << 
            number_subregions << ": " << std::endl;

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

    std::vector<std::vector<int>> getCandidatePatterns(const std::vector<std::vector<int>>& prevalentPattern, int degree, const std::map<int, FeatureInformation> featureInfo) {
        std::set<std::vector<int>> prevalentPatterns(prevalentPattern.begin(), prevalentPattern.end());
        // Extract features from the keys of featureInfo
        std::vector<int> features;
        for (const auto& pair : featureInfo) {
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

    std::vector<std::vector<int>> degree2Processing(std::vector<std::vector<int>> candidatePatterns, 
                                                    double prevalence_threshold, int id, 
                                                    std::map<int, FeatureInformation> featureInfo,
                                                    std::map<int, std::vector<int>> starNeighbors,
                                                    std::map<std::vector<int>, std::map<std::vector<int>,
                                                        std::vector<int>>> &instance_table, 
                                                    std::map<std::vector<int>,
                                                        std::map<int, std::set<int>>> &hashmap) {   
        std::vector<std::vector<int>> size2_patterns;

        for (const auto& coloc : candidatePatterns) {
            int first_feature = coloc[0];
            int second_feature = coloc[1];

            // determine first feature start and end
            auto feature1 = featureInfo.find(first_feature);
            FeatureInformation values1 = feature1->second;
            int first_feature_start = values1.start;
            int first_feature_end = values1.end;
            // determine second feature start and end
            auto feature2 = featureInfo.find(second_feature);
            FeatureInformation values2 = feature2->second;
            int second_feature_start = values2.start;
            int second_feature_end = values2.end;

            std::vector<int> coloc_key = {first_feature, second_feature};
            instance_table[coloc_key] = {};
            hashmap[coloc_key] = {};
            hashmap[coloc_key][first_feature] = {};
            hashmap[coloc_key][second_feature] = {};

            for (int index = first_feature_start; index <= first_feature_end; index++) {
                auto star_neighbor_it = starNeighbors.find(index);
                std::vector<int> neighbors = findNeighborsInRange(star_neighbor_it->second, 
                                                                  second_feature_start,
                                                                  second_feature_end);
                if (!neighbors.empty()) {
                    std::vector<int> index_tuple = {index};
                    instance_table[coloc_key][index_tuple] = neighbors;
                    hashmap[coloc_key][first_feature].insert(index);
                    for (int neighbor : neighbors) {
                        hashmap[coloc_key][second_feature].insert(neighbor);
                    }
                }
            }

            double pr_first_feature = static_cast<double>(hashmap[coloc_key][first_feature].size()) /
                featureInfo[first_feature].count;
            double pr_second_feature = static_cast<double>(hashmap[coloc_key][second_feature].size()) /
                featureInfo[second_feature].count;

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
        }

        std::cout << "Degree 2 Prevalent Patterns for Sub-Region " << id << ":" << std::endl;
        for (auto i : size2_patterns) {
            std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
        }
        return size2_patterns;
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
    
    
    
    void subregion_main() {
        std::string input_file_path = "/home/amk7r/colocation_mining/updated_regional/Data/Region6/";
        std::string distance_file = 
         "/home/amk7r/colocation_mining/updated_regional/IntermediateData/distance_threshold.txt";
        //distance_threshold = read_distance(distance_file);
        std::vector<int> offsets;
        offsets.push_back(0);
        std::vector<int> indices;
        std::vector<Feature> features = read_data(input_file_path, offsets, indices);
        offsets.pop_back();
        
        int number_subregions = offsets.size();
        for(int i = 0; i < number_subregions; i++) {
            subregions.push_back(Subregion(i));
        }
        
        int k = 0;
        for(int i = 0; i < number_subregions; i++) {
            subregions[i].data.assign(features.begin() + offsets[i] + indices[k], 
                                  features.begin() + offsets[i] + indices[k + 1] + 1);        
            subregions[i].featureInfo = calc_featureInfo(features, indices[k], 
                                                         indices[k + 1], offsets[i]);
            k += 2;

            subregions[i].starNeighbors = calc_starNeighbors(subregions[i].data);

            // find prevalent patterns in subregions
            std::vector<std::vector<int>> size2_candidatePatterns = generate_size2_combos(
                                                                subregions[i].featureInfo);
            subregions[i].size2_patterns = degree2Processing(size2_candidatePatterns,
                                                             prevalence_threshold, i, 
                                                             subregions[i].featureInfo,
                                                             subregions[i].starNeighbors,
                                                             subregions[i].instance_table,
                                                             subregions[i].hashmap);
            int degree = 3;
            std::vector<std::vector<int>> candidatePatterns =
                getCandidatePatterns(subregions[i].size2_patterns, degree,
                                     subregions[i].featureInfo);
            while (!candidatePatterns.empty()) {
                std::vector<std::vector<int>> prevalent_patterns = 
                colocationGeneral(candidatePatterns, prevalence_threshold, degree, i,
                                 subregions[i].instance_table, subregions[i].hashmap,
                                 subregions[i].starNeighbors, subregions[i].featureInfo);
                degree += 1;
                if (prevalent_patterns.size() == 0) {
                    break;
                }
                
                if(degree == 5) {
                    break;
                }
                candidatePatterns = getCandidatePatterns(prevalent_patterns,
                                                         degree, subregions[i].featureInfo);
            }
        }  
    }
    
    std::vector<Feature> read_border_data(const std::string& folderPath) {
        std::vector<Feature> features;
        int region_count = 0;
        int row_count = 0;
        int total_count = 0;


        std::ifstream file(folderPath);
        if (!file.is_open()) {
            std::cerr << "Error opening border file!" << std::endl;
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
                else if (colIndex == 3) {
                    tempFeature.id = std::stoi(value);
                }
                colIndex++;
            }
            
            features.push_back(tempFeature);
        }
        file.close();

        return features;
    }
    
    std::map<int, FeatureInformation> border_calc_featureInfo(std::vector<Feature> features) {
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
    
    std::map<int, std::vector<int>> border_calc_starNeighbors(std::vector<Feature> features) {
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

    void border_main(int number_borders) {
        for (int i = 0; i < number_borders; i++) {
            borders.push_back(Border(i));
        }
        
        std::string input_file = "/home/amk7r/colocation_mining/updated_regional/IntermediateData/border.csv";
        std::vector<Feature> features = read_border_data(input_file);
        
        for (int i = 0; i < number_borders; i++) {
            borders[i].featureInfo = border_calc_featureInfo(features);
            borders[i].starNeighbors = border_calc_starNeighbors(features);
            std::vector<std::vector<int>> size2_candidatePatterns = generate_size2_combos(
                borders[i].featureInfo);
            borders[i].size2_patterns = degree2Processing(size2_candidatePatterns, prevalence_threshold, i, 
                                                          borders[i].featureInfo, borders[i].starNeighbors,
                                                          borders[i].instance_table, borders[i].hashmap);

        } 
        return;
    }
    
    void update_border_info(int* ids, int ids_size, int i) {
        // update the border hashmap so that it has the original IDS not the indicies
        std::map<std::vector<int>, std::map<int, std::set<int>>> temp_hash;
        for (const auto& outer_pair : borders[i].hashmap) {
            const std::vector<int>& outer_key = outer_pair.first;
            const auto& inner_map = outer_pair.second;
            std::map<int, std::set<int>> temp_inner_map;
            for (const auto& inner_pair : inner_map) {
                const int& inner_key = inner_pair.first;
                const std::set<int>& inner_set = inner_pair.second;
                std::set<int> temp_inner_set;
                for (int index : inner_set) {
                    if (index < ids_size) {
                        temp_inner_set.insert(ids[index]);
                    } else {
                        std::cout << "Index out of range: " << index << std::endl;
                    }
                }
                temp_inner_map[inner_key] = temp_inner_set;
            }
            temp_hash[outer_key] = temp_inner_map;
        }
        borders[i].hashmap = temp_hash;
    
        // update the border instance table so that it has the original IDS not indicies
        std::map<std::vector<int>, std::map<std::vector<int>, std::vector<int>>> temp_instance_table;

        for (const auto& outer_pair : borders[i].instance_table) {
            const std::vector<int>& outer_key = outer_pair.first;
            const auto& inner_map = outer_pair.second;
            std::map<std::vector<int>, std::vector<int>> temp_inner_map;
            for (const auto& inner_pair : inner_map) {
                const std::vector<int>& inner_key = inner_pair.first;
                const std::vector<int>& inner_value = inner_pair.second;
                std::vector<int> new_inner_key;
                for (int index : inner_key) {
                    if (index < ids_size) {
                        new_inner_key.push_back(ids[index]);
                    } else {
                        std::cout << "Index out of range: " << index << std::endl;
                    }
                }
                std::vector<int> new_inner_value;
                for (int index : inner_value) {
                    if (index < ids_size) {
                        new_inner_value.push_back(ids[index]);
                    } else {
                        std::cout << "Index out of range: " << index << std::endl;
                    }
                }
                temp_inner_map[new_inner_key] = new_inner_value;
            }
            temp_instance_table[outer_key] = temp_inner_map;
        }
        borders[i].instance_table = temp_instance_table;
        
        // update the border star neighbors so that is has the original IDS not indices
        std::map<int, std::vector<int>> temp_star_neighbors;
        
        for (const auto& entry : borders[i].starNeighbors) {
            int original_id = ids[entry.first]; // Get the original ID corresponding to the index
            const std::vector<int>& neighbors_indices = entry.second;
            std::vector<int> original_neighbors;
            for (int index : neighbors_indices) {
                if (index < ids_size) {
                    original_neighbors.push_back(ids[index]); // Get the original ID for each neighbor
                } else {
                    std::cout << "Index out of range: " << index << std::endl;
                }
            }
            temp_star_neighbors[original_id] = original_neighbors;
        }
        borders[i].starNeighbors = temp_star_neighbors;        
    }
    
    void combine_hashmaps(int number_subregions, int number_borders) {
        std::set<std::vector<int>> keys_needed;
        
        // determine which keys from the subregions and borders need to be combined to form the regional hashmap
        for (int i = 0; i < number_subregions; i++) {
            std::set<std::vector<int>> temp_vect(subregions[i].size2_patterns.begin(), 
                                                         subregions[i].size2_patterns.end());  
            if (i == 0) {
                keys_needed = temp_vect;
            } else {
                std::set<std::vector<int>> temp_result;
                std::set_union(keys_needed.begin(), keys_needed.end(),
                               temp_vect.begin(), temp_vect.end(),
                               std::inserter(temp_result, temp_result.begin()));
                keys_needed = temp_result;
            }
        }
        
        for (int i = 0; i < number_borders; i++) {
            std::set<std::vector<int>> temp_vect(borders[i].size2_patterns.begin(), 
                                                         borders[i].size2_patterns.end());
            std::set<std::vector<int>> temp_result;
            std::set_union(keys_needed.begin(), keys_needed.end(),
                           temp_vect.begin(), temp_vect.end(),
                           std::inserter(temp_result, temp_result.begin()));
            keys_needed = temp_result;
        }
        
        // combine the hashmaps
        for (const auto& key : keys_needed) {
            for (int i = 0; i < number_subregions; ++i) {
                auto it = subregions[i].hashmap.find(key);
                if (it != subregions[i].hashmap.end()) {
                    for (const auto& inner_pair : it->second) {
                        const auto& inner_key = inner_pair.first;
                        const auto& inner_value = inner_pair.second;
                        auto& region_inner_map = region.hashmap[key][inner_key];
                        region_inner_map.insert(inner_value.begin(), inner_value.end());
                    }
                }
            }

            for (int i = 0; i < number_borders; ++i) {
                auto it = borders[i].hashmap.find(key);
                if (it != borders[i].hashmap.end()) {
                    for (const auto& inner_pair : it->second) {
                        const auto& inner_key = inner_pair.first;
                        const auto& inner_value = inner_pair.second;
                        auto& region_inner_map = region.hashmap[key][inner_key];
                        region_inner_map.insert(inner_value.begin(), inner_value.end());
                    }
                }
            }
        }  
    }
    
    void combine_instance_tables(int number_subregions, int number_borders) {
        std::set<std::vector<int>> keys_needed;
        
        // determine which keys from the subregions and borders need to be combined to form the regional hashmap
        for (int i = 0; i < number_subregions; i++) {
            std::set<std::vector<int>> temp_vect(subregions[i].size2_patterns.begin(), 
                                                         subregions[i].size2_patterns.end());  
            if (i == 0) {
                keys_needed = temp_vect;
            } else {
                std::set<std::vector<int>> temp_result;
                std::set_union(keys_needed.begin(), keys_needed.end(),
                               temp_vect.begin(), temp_vect.end(),
                               std::inserter(temp_result, temp_result.begin()));
                keys_needed = temp_result;
            }
        }
        
        for (int i = 0; i < number_borders; i++) {
            std::set<std::vector<int>> temp_vect(borders[i].size2_patterns.begin(), 
                                                         borders[i].size2_patterns.end());
            std::set<std::vector<int>> temp_result;
            std::set_union(keys_needed.begin(), keys_needed.end(),
                           temp_vect.begin(), temp_vect.end(),
                           std::inserter(temp_result, temp_result.begin()));
            keys_needed = temp_result;
        }
        
        // Combine the instance tables
        for (const auto& key : keys_needed) {
            for (int i = 0; i < number_subregions; ++i) {
                auto it = subregions[i].instance_table.find(key);
                if (it != subregions[i].instance_table.end()) {
                    for (const auto& inner_pair : it->second) {
                        const auto& inner_key = inner_pair.first;
                        const auto& inner_value = inner_pair.second;
                        auto& region_inner_map = region.instance_table[key];
                        if (region_inner_map.find(inner_key) != region_inner_map.end()) {
                            region_inner_map[inner_key].insert(region_inner_map[inner_key].end(), inner_value.begin(),
                                                               inner_value.end());
                        } else {
                            region_inner_map[inner_key] = inner_value;
                        }
                    }
                }
            }

            for (int i = 0; i < number_borders; ++i) {
                auto it = borders[i].instance_table.find(key);
                if (it != borders[i].instance_table.end()) {
                    for (const auto& inner_pair : it->second) {
                        const auto& inner_key = inner_pair.first;
                        const auto& inner_value = inner_pair.second;
                        auto& region_inner_map = region.instance_table[key];
                        if (region_inner_map.find(inner_key) != region_inner_map.end()) {
                            region_inner_map[inner_key].insert(region_inner_map[inner_key].end(), inner_value.begin(),
                                                               inner_value.end());
                        } else {
                            region_inner_map[inner_key] = inner_value;
                        }
                    }
                }
            }
        }
    }
    
    std::vector<std::vector<int>> region_generate_size2_combos() {
        std::set<std::vector<int>> patterns;
        for (const auto& entry : region.hashmap) {
            patterns.insert(entry.first);
        }
        std::vector<std::vector<int>> size2_candidatePatterns(patterns.begin(), patterns.end());
        return size2_candidatePatterns;
    }
    
    std::vector<std::vector<int>> region_degree2Processing(std::vector<std::vector<int>> candidatePatterns, 
                                                                int candidatePatterns_size, 
                                                            double prevalence_threshold,
                                                                int number_subregions) {
            std::vector<std::vector<int>> prevalent;
            
            for (const auto& coloc : candidatePatterns) {
                int first_feature = coloc[0];
                int second_feature = coloc[1];
                int total_first_feature = 0;
                int total_second_feature = 0;
                
                for (int i = 0; i < number_subregions; i++) {
                    auto feature1 = subregions[i].featureInfo.find(first_feature);
                    FeatureInformation values1 = feature1->second;
                    total_first_feature += values1.count;
                    
                    auto feature2 = subregions[i].featureInfo.find(second_feature);
                    FeatureInformation values2 = feature2->second;
                    total_second_feature += values2.count;
                }
                
                double pr_first_feature = static_cast<double>(region.hashmap[coloc][first_feature].size()) /
                    total_first_feature;
                double pr_second_feature = static_cast<double>(region.hashmap[coloc][second_feature].size()) /
                    total_second_feature;
                //std::cout << this->hashmap[coloc][first_feature].size() << " / " << total_first_feature << std::endl;
                //std::cout << this->hashmap[coloc][second_feature].size() << " / " << total_second_feature << std::endl;
                double PI = 0.0;

                if (pr_first_feature < pr_second_feature) {
                    PI = pr_first_feature;
                } else {
                    PI = pr_second_feature;
                }

                if (PI >= prevalence_threshold) {
                       prevalent.push_back(coloc);
                }
            }
            
            std::cout << "Degree 2 Prevalent Patterns for Entire Region:" << std::endl;
            for (auto i : prevalent) {
                std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
            }
            return prevalent;
        }
    
    std::vector<std::vector<int>> region_getCandidatePatterns(
        const std::vector<std::vector<int>>& prevalentPattern, int degree, 
        const std::vector<int> features) {
        std::set<std::vector<int>> prevalentPatterns(prevalentPattern.begin(), prevalentPattern.end());

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
    
    std::vector<std::vector<int>> region_colocationGeneral(
        std::vector<std::vector<int>> candidatePatterns, int candidatePatterns_size, 
        double prevalence_threshold, int degree, int number_subregions) {
        
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
                    std::set<int> lastFeatureNeighbors;
                    
                    /*for (auto subregion : subregions) {
                        int lastFeatureStart = subregion.featureInfo[lastFeature].start;
                        int lastFeatureEnd = subregion.featureInfo[lastFeature].end;

                        auto subregion_star_neighbor_it =
                            subregion.starNeighbors.find(instanceID);
                        if (subregion_star_neighbor_it != subregion.starNeighbors.end()) {
                            std::vector<int> subregion_star_neighbor_list(
                                subregion_star_neighbor_it->second.begin(), 
                                subregion_star_neighbor_it->second.end());
                            std::vector<int> neighborsInRange_vector = 
                                findNeighborsInRange(subregion_star_neighbor_list,
                                                     lastFeatureStart, 
                                                     lastFeatureEnd);
                            std::set<int> neighborsInRange(neighborsInRange_vector.begin(),
                                                           neighborsInRange_vector.end());
                            lastFeatureNeighbors.insert(neighborsInRange.begin(),
                                                        neighborsInRange.end());
                        }
                        auto border_star_neighbor_it = 
                            borders[0].starNeighbors.find(instanceID);
                        if (border_star_neighbor_it != borders[0].starNeighbors.end()) {
                            std::vector<int> border_star_neighbor_list(
                                border_star_neighbor_it->second.begin(), 
                                border_star_neighbor_it->second.end());
                            std::vector<int> neighborsInRange_vector = 
                                findNeighborsInRange(border_star_neighbor_list,
                                                     lastFeatureStart, 
                                                     lastFeatureEnd);
                            std::set<int> neighborsInRange(neighborsInRange_vector.begin(),
                                                           neighborsInRange_vector.end());
                            lastFeatureNeighbors.insert(neighborsInRange.begin(),
                                                        neighborsInRange.end());
                        }
                    }*/
                    
                    for (const auto& subregion : subregions) {
                        const int lastFeatureStart = subregion.featureInfo.at(lastFeature).start;
                        const int lastFeatureEnd = subregion.featureInfo.at(lastFeature).end;

                        // Define a lambda function to process neighbors
                        auto processNeighbors = [&](const auto& neighborMap) {
                            auto neighbor_it = neighborMap.find(instanceID);
                            if (neighbor_it != neighborMap.end()) {
                                const auto& neighbor_list = neighbor_it->second; // Use reference
                                auto neighborsInRange_vector = findNeighborsInRange(neighbor_list, lastFeatureStart, lastFeatureEnd);
                                lastFeatureNeighbors.insert(neighborsInRange_vector.begin(), neighborsInRange_vector.end());
                            }
                        };

                        // Process subregion neighbors
                        processNeighbors(subregion.starNeighbors);

                        // Process border neighbors
                        processNeighbors(borders[0].starNeighbors);
                    }
                    
                    if (commonLastNeighbors.empty()) {
                        commonLastNeighbors = lastFeatureNeighbors;
                    } else {
                        std::set<int> intersection;
                        std::set_intersection(commonLastNeighbors.begin(),
                                              commonLastNeighbors.end(),
                                              lastFeatureNeighbors.begin(),
                                              lastFeatureNeighbors.end(),
                                              std::inserter(intersection,
                                              intersection.begin()));
                        commonLastNeighbors = intersection;
                    }
                }
                for (int n : colocTable[key]) {
                    std::set<int> all_neighbors;
                    /*for (const auto& subregion : subregions) {
                        const int lastFeatureStart = subregion.featureInfo.at(lastFeature).start;
                        const int lastFeatureEnd = subregion.featureInfo.at(lastFeature).end;
                        
                        std::set<int> neighbors_subregion;
                        
                        auto subregion_star_neighbor_it = subregion.starNeighbors.find(n);
                        if (subregion_star_neighbor_it != subregion.starNeighbors.end()) {
                            std::vector<int> subregion_star_neighbor_list(
                                subregion_star_neighbor_it->second.begin(), 
                                subregion_star_neighbor_it->second.end());
                            std::vector<int> neighborsInRange_vector =
                                findNeighborsInRange(subregion_star_neighbor_list,
                                                     lastFeatureStart, 
                                                     lastFeatureEnd);
                            std::set<int> neighborsInRange(neighborsInRange_vector.begin(),
                                                           neighborsInRange_vector.end());
                            std::set_intersection(commonLastNeighbors.begin(),
                                                  commonLastNeighbors.end(),
                                                  neighborsInRange.begin(),
                                                  neighborsInRange.end(),
                                                  std::inserter(neighbors_subregion, 
                                                  neighbors_subregion.begin()));
                        }
                        std::set<int> neighbors_border;
                        auto border_star_neighbor_it = borders[0].starNeighbors.find(n);
                        if (border_star_neighbor_it != borders[0].starNeighbors.end()) {
                            std::vector<int> border_star_neighbor_list(
                                border_star_neighbor_it->second.begin(), 
                                border_star_neighbor_it->second.end());
                            std::vector<int> neighborsInRange_vector = 
                                findNeighborsInRange(border_star_neighbor_list,
                                                     lastFeatureStart, 
                                                     lastFeatureEnd);
                            std::set<int> neighborsInRange(neighborsInRange_vector.begin(),
                                                           neighborsInRange_vector.end());
                                                           std::set_intersection(
                                                           commonLastNeighbors.begin(),
                                                           commonLastNeighbors.end(),
                                                           neighborsInRange.begin(), 
                                                           neighborsInRange.end(),
                                                           std::inserter(neighbors_border,
                                                           neighbors_border.begin()));
                        }
                        all_neighbors.insert(neighbors_subregion.begin(),
                                             neighbors_subregion.end());
                        all_neighbors.insert(neighbors_border.begin(), 
                                             neighbors_border.end());
                    }*/
                    
                    
                    for (const auto& subregion : subregions) {
                        const int lastFeatureStart = subregion.featureInfo.at(lastFeature).start;
                        const int lastFeatureEnd = subregion.featureInfo.at(lastFeature).end;

                        // Use a single set to collect neighbors from both sources
                        std::set<int> neighbors_subregion;

                        // Process subregion star neighbors
                        auto subregion_star_neighbor_it = subregion.starNeighbors.find(n);
                        if (subregion_star_neighbor_it != subregion.starNeighbors.end()) {
                            const auto& subregion_star_neighbor_list = subregion_star_neighbor_it->second; // Use reference
                            auto neighborsInRange_vector = findNeighborsInRange(subregion_star_neighbor_list, lastFeatureStart, lastFeatureEnd);
                            std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());

                            // Compute intersection directly into neighbors_subregion
                            std::set_intersection(
                                commonLastNeighbors.begin(),
                                commonLastNeighbors.end(),
                                neighborsInRange.begin(),
                                neighborsInRange.end(),
                                std::inserter(neighbors_subregion, neighbors_subregion.begin())
                            );
                        }

                        // Process border star neighbors
                        std::set<int> neighbors_border; // Only need to declare this once
                        auto border_star_neighbor_it = borders[0].starNeighbors.find(n);
                        if (border_star_neighbor_it != borders[0].starNeighbors.end()) {
                            const auto& border_star_neighbor_list = border_star_neighbor_it->second; // Use reference
                            auto neighborsInRange_vector = findNeighborsInRange(border_star_neighbor_list, lastFeatureStart, lastFeatureEnd);
                            std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());

                            // Compute intersection directly into neighbors_border
                            std::set_intersection(
                                commonLastNeighbors.begin(),
                                commonLastNeighbors.end(),
                                neighborsInRange.begin(),
                                neighborsInRange.end(),
                                std::inserter(neighbors_border, neighbors_border.begin())
                            );
                        }

                        // Insert the results into all_neighbors set
                        all_neighbors.insert(neighbors_subregion.begin(), neighbors_subregion.end());
                        all_neighbors.insert(neighbors_border.begin(), neighbors_border.end());
                    }
                    
                    if (!all_neighbors.empty()) {
                        std::vector<int> new_key = key;
                        new_key.push_back(n);
                        std::vector<int> intersectionVec(all_neighbors.begin(), 
                                                         all_neighbors.end());
                        region.instance_table[currentPattern][new_key] = intersectionVec;

                        for (size_t k = 0; k < new_key.size(); ++k) {
                            region.hashmap[currentPattern][currentPattern[k]].insert(
                                new_key[k]);
                        }
                        region.hashmap[currentPattern][lastFeature].insert(
                            all_neighbors.begin(), all_neighbors.end());
                    }
                }
            }
            std::vector<double> pr;
                for (int m = 0; m < degree; ++m) {
                    int f = currentPattern[m];
                    int total_count = 0;
                    for (auto subregion : subregions) {
                        total_count += subregion.featureInfo[f].count;
                    }
                    double ratio = static_cast<double>(
                        region.hashmap[currentPattern][f].size()) / total_count;
                    std::cout << region.hashmap[currentPattern][f].size() << " / " << total_count << std::endl;
                    pr.push_back(ratio);
                }
                double PI = *std::min_element(pr.begin(), pr.end());
                if (PI >= prevalence_threshold) {
                    prevalent.push_back(currentPattern);
                }
        }
        std::cout << "Degree " << degree << " Prevalent Patterns for Entire Region:" +
            std::to_string(number_subregions) + ":"<< std::endl;
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
        return prevalent;
    }
    
    void region_main() {
        std::set<int> set_features;
        for(int i = 0; i < subregions.size(); i++) {
            for(const auto &[key, _] : subregions[i].featureInfo) {
                set_features.insert(key);
            }
        }
        
        for(int i = 0; i < borders.size(); i++) {
            for(const auto &[key, _] : borders[i].featureInfo) {
                set_features.insert(key);
            }
        }
                
        std::vector<std::vector<int>> size2_candidatePatterns = 
            region_generate_size2_combos();
        std::vector<std::vector<int>> prevalentPatterns = 
            region_degree2Processing(size2_candidatePatterns,                               
                                     size2_candidatePatterns.size(),                       
                                     prevalence_threshold, subregions.size());
        
        int degree = 3;
        std::vector<int> features(set_features.begin(), set_features.end());
        std::vector<std::vector<int>> candidatePatterns =
            region_getCandidatePatterns(prevalentPatterns, degree, features);
        while (!candidatePatterns.empty()) {
            std::vector<std::vector<int>> prevalentPatterns =
                region_colocationGeneral(candidatePatterns, candidatePatterns.size(), 
                                         prevalence_threshold, degree, subregions.size());
            degree += 1;
            if (prevalentPatterns.size() == 0) {
                break;
            }
            
            if (degree == 5) {
                break;
            }
            
            candidatePatterns = region_getCandidatePatterns(prevalentPatterns, degree,
                                                            features);
        }
        
    }

}