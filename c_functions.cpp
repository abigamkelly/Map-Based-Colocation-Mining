#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>

extern "C" {
    // Holds all information pertaining to the sub-regions
    class SubRegion {
    public:
        int subregion_id;  // unique identifier
        std::map<std::string, std::vector<int>> featureInfo;
        std::map<int, std::vector<int>> star_neighbors;
        std::vector<std::vector<std::string>> size2_patterns;
        std::map<std::vector<std::string>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<std::string>, std::map<std::string, std::set<int>>> hashmap;
        
        SubRegion(int number) {
            subregion_id = number;
        }

        // get the featureInfo
        void read_featureInfo() {
            std::ifstream featureInfo_file("required_files/featureInfo/featureInfo" + std::to_string(this->subregion_id) + ".csv");
        
            if (!featureInfo_file.is_open()) {
                std::cerr << "Error opening file!" << std::endl;
            }

            std::string line;
            std::getline(featureInfo_file, line);

            while (std::getline(featureInfo_file, line)) {
                std::istringstream lineStream(line);
                std::string cell;
                std::string key;
                std::vector<int> values;
                if (std::getline(lineStream, cell, ',')) {
                    key = cell;
                }
                while (std::getline(lineStream, cell, ',')) {
                    values.push_back(std::stoi(cell));
                }
                this->featureInfo[key] = values;
            }

            featureInfo_file.close();
        }

        // get the star_neighbors
        void read_star_neighbors() {
            std::ifstream star_neighbors_file("required_files/starNeighbors/starNeighbors" + 
                                              std::to_string(this->subregion_id) + ".csv");

            if (!star_neighbors_file.is_open()) {
                std::cerr << "Error opening file" << std::endl;
            }
            std::string line;
            std::getline(star_neighbors_file, line);

            while (std::getline(star_neighbors_file, line)) {
                std::stringstream ss(line);
                int key;
                ss >> std::ws;
                ss >> key;
                ss.ignore();

                int value;
                std::vector<int> values;
                while (ss >> value) {
                    values.push_back(value);
                    ss.ignore();
                }
                star_neighbors[key] = values;
            }
            star_neighbors_file.close();
        }

        // generate the size 2 candidate patterns
        std::vector<std::vector<std::string>> generate_size2_combos() {
            std::vector<std::string> features;
            for (const auto &entry : this->featureInfo) {
                features.push_back(entry.first);
            }

            std::vector<std::vector<std::string>> size2_candidatePatterns;
            for (size_t i = 0; i < features.size(); ++i) {
                for (size_t j = i + 1; j < features.size(); ++j) {
                    size2_candidatePatterns.push_back({features[i], features[j]});
                }
            }
            return size2_candidatePatterns;
        }

        // log(n) search that finds neighbors within a specific start and end range
        std::vector<int> findNeighborsInRange(const std::vector<int>& arr, int x, int y) {
            auto start_it = std::lower_bound(arr.begin(), arr.end(), x);
            auto end_it = std::upper_bound(arr.begin(), arr.end(), y);
            return std::vector<int>(start_it, end_it);
        }

        // find the degree 2 prevalent patterns
        std::vector<std::vector<std::string>> degree2Processing(
            std::vector<std::vector<std::string>> candidatePatterns, int candidatePatterns_size, 
            double prevalence_threshold, int number_subregions) {   
            std::vector<std::vector<std::string>> size2_patterns;
            
            for (const auto& coloc : candidatePatterns) {
                std::string first_feature = coloc[0];
                std::string second_feature = coloc[1];
                
                // determine first feature start and end
                auto feature1 = this->featureInfo.find(first_feature);
                std::vector<int> values1 = feature1->second;
                int first_feature_start = values1[1];
                int first_feature_end = values1[2];
                // determine second feature start and end
                auto feature2 = this->featureInfo.find(second_feature);
                std::vector<int> values2 = feature2->second;
                int second_feature_start = values2[1];
                int second_feature_end = values2[2];
                
                std::vector<std::string> coloc_key = {first_feature, second_feature};
                instance_table[coloc_key] = {};
                hashmap[coloc_key] = {};
                hashmap[coloc_key][first_feature] = {};
                hashmap[coloc_key][second_feature] = {};
                
                for (int index = first_feature_start; index <= first_feature_end; index++) {
                    auto star_neighbor_it = star_neighbors.find(index);
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

                // calculate the participation rations
                double pr_first_feature = static_cast<double>(hashmap[coloc_key][first_feature].size()) /
                    featureInfo[first_feature][0];
                double pr_second_feature = static_cast<double>(hashmap[coloc_key][second_feature].size()) /
                    featureInfo[second_feature][0];
                double PI = 0.0;

                // calculate the participation index
                if (pr_first_feature < pr_second_feature) {
                    PI = pr_first_feature;
                } else {
                    PI = pr_second_feature;
                }

                if (PI >= prevalence_threshold) {
                       size2_patterns.push_back(coloc_key);
                }
            }
            
            std::cout << "Degree 2 Prevalent Patterns for Sub-Region " + std::to_string(number_subregions) 
                + ":" << std::endl;
            std::ofstream size2_file("size2_patterns_subregion" + std::to_string(number_subregions) + ".txt");
            for (auto i : size2_patterns) {
                std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
                size2_file << "(" << i[0] << ", " << i[1] << ")" << std::endl;
            }
            size2_file.close();
            return size2_patterns;
        }

        // generates the combinations for the candidate patterns with degree > 2
        void generateCombinations(const std::vector<std::string>& features, int degree, 
                                  std::vector<std::vector<std::string>>& result, 
                                  std::vector<std::string>& current, int start) {
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

        // check if all (degree-1)-subpatterns of a pattern are in the prevalent patterns
        bool allSubpatternsInPrevalent(const std::vector<std::string>& pattern, 
                                       const std::set<std::vector<std::string>>& prevalentPatterns, 
                                       int degree) {
            std::vector<std::vector<std::string>> subpatterns;
            std::vector<std::string> current;
            generateCombinations(pattern, degree - 1, subpatterns, current, 0);

            for (const auto& subpattern : subpatterns) {
                if (prevalentPatterns.find(subpattern) == prevalentPatterns.end()) {
                    return false;
                }
            }
            return true;
        }

        std::vector<std::vector<std::string>> getCandidatePatterns(
            const std::vector<std::vector<std::string>>& prevalentPattern, int degree) {
            std::set<std::vector<std::string>> prevalentPatterns(prevalentPattern.begin(), prevalentPattern.end());
            // Extract features from the keys of featureInfo
            std::vector<std::string> features;
            for (const auto& pair : this->featureInfo) {
                features.push_back(pair.first);
            }

            std::vector<std::vector<std::string>> _candidatePatterns;
            std::vector<std::vector<std::string>> _patterns;
            std::vector<std::string> current;

            generateCombinations(features, degree, _patterns, current, 0);

            for (const auto& pattern : _patterns) {
                if (allSubpatternsInPrevalent(pattern, prevalentPatterns, degree)) {
                    _candidatePatterns.push_back(pattern);
                }
            }

            return _candidatePatterns;
        }

        // find prevalent patterns for degree > 2
        std::vector<std::vector<std::string>> colocationGeneral(
            std::vector<std::vector<std::string>> candidatePatterns, int candidatePatterns_size, 
                           double prevalence_threshold, int degree, int number_subregions) {
            
            std::vector<std::vector<std::string>> prevalent;
            
            for (const auto& currentPattern : candidatePatterns) {
                std::vector<std::string> basePattern;
                for (int j = 0; j < degree - 1; j++) {
                    basePattern.push_back(currentPattern[j]);
                }
                std::string lastFeature = currentPattern[degree - 1];
                // Add entries to instance_table and hashmap
                this->instance_table[currentPattern] = {};
                this->hashmap[currentPattern] = {};
                
                // Initialize sets for each element in currentPattern in hashmap
                for (const auto& f : currentPattern) {
                    this->hashmap[currentPattern][f] = {};
                }
                
                auto colocTableIt = this->instance_table.find(basePattern);
                std::map<std::vector<int>, std::vector<int>>& colocTable = colocTableIt->second;
                
                for (const auto& entry : colocTable) {
                    const std::vector<int>& key = entry.first;
                    std::set<int> commonLastNeighbors;
                    for (int instanceID : key) {
                        auto star_neighbor_it = this->star_neighbors.find(instanceID);
                        if (commonLastNeighbors.empty()) {
                            std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second,
                                                                         this->featureInfo[lastFeature][1],
                                                                         this->featureInfo[lastFeature][2]);
                            commonLastNeighbors.insert(temp_vector.begin(), temp_vector.end());
                        } else {
                            std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second, 
                                                                                  this->featureInfo[lastFeature][1],
                                                                                  this->featureInfo[lastFeature][2]);
                            std::set<int> temp_commonLastNeighbors(temp_vector.begin(), temp_vector.end());
                            std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                                  temp_commonLastNeighbors.begin(), temp_commonLastNeighbors.end(),
                                                  std::inserter(commonLastNeighbors, commonLastNeighbors.begin()));
                        }
                    }
                    for (int n : colocTable[key]) {
                        auto star_neighbor_it = this->star_neighbors.find(n);
                        std::vector<int> temp_vect = findNeighborsInRange(star_neighbor_it->second,
                                                                          this->featureInfo[lastFeature][1],
                                                                          this->featureInfo[lastFeature][2]);
                        std::set<int> temp_neighbors(temp_vect.begin(), temp_vect.end());
                        std::set<int> neighbors;
                        std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                              temp_neighbors.begin(), temp_neighbors.end(),
                                              std::inserter(neighbors, neighbors.begin()));

                        if (!neighbors.empty()) {
                            std::vector<int> new_key = key;
                            new_key.push_back(n);
                            std::vector<int> intersectionVec(neighbors.begin(), neighbors.end());
                            this->instance_table[currentPattern][new_key] = intersectionVec;

                            for (size_t k = 0; k < new_key.size(); ++k) {
                                this->hashmap[currentPattern][currentPattern[k]].insert(new_key[k]);
                            }
                            this->hashmap[currentPattern][lastFeature].insert(neighbors.begin(), neighbors.end());
                        }
                    }
                }
            
                std::vector<double> pr;
                for (int m = 0; m < degree; ++m) {
                    std::string f = currentPattern[m];
                    double ratio = static_cast<double>(this->hashmap[currentPattern][f].size()) 
                        / this->featureInfo[f][0];
                    pr.push_back(ratio);
                }
                double PI = *std::min_element(pr.begin(), pr.end());
                if (PI >= prevalence_threshold) {
                    prevalent.push_back(currentPattern);
                }
            }

            std::cout << "Degree " << degree << " Prevalent Patterns for Sub-Region:" +
                std::to_string(number_subregions) + ":"<< std::endl;
            std::ofstream patterns_file("size" + std::to_string(degree) + "_patterns_subregion" 
                                        + std::to_string(number_subregions) + ".txt");
            
            for (auto& patternVec : prevalent) {
                std::cout << "(";
                patterns_file << "(";
                for (size_t i = 0; i < patternVec.size(); i++) {
                    std::cout << patternVec[i];
                    patterns_file << patternVec[i];
                    if (i < patternVec.size() - 1) {
                        std::cout << ", ";
                        patterns_file << ", ";
                    }
                }
                std::cout << ")" << std::endl;
                patterns_file << ")" << std::endl;
            }
            patterns_file.close();
            return prevalent;
        }
    };
    std::vector<SubRegion> subregions;

    // this class holds all information pertaining to the border region
    class Border {
    public:
        int border_id;  // unique identifier
        std::map<std::string, std::vector<int>> featureInfo;
        std::map<int, std::vector<int>> star_neighbors;
        std::vector<std::vector<std::string>> size2_patterns;
        std::map<std::vector<std::string>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<std::string>, std::map<std::string, std::set<int>>> hashmap;
        
        Border(int number) {
            border_id = number;
        }
        
        void read_featureInfo() {
            std::ifstream featureInfo_file("required_files/border_featureInfo/featureInfo.csv");
        
            if (!featureInfo_file.is_open()) {
                std::cerr << "Error opening file" << std::endl;
            }

            std::string line;
            std::getline(featureInfo_file, line);  // Skip the first row (header)

            while (std::getline(featureInfo_file, line)) {
                std::istringstream lineStream(line);
                std::string cell;
                std::string key;
                std::vector<int> values;

                // Read the first cell (key)
                if (std::getline(lineStream, cell, ',')) {
                    key = cell; // Assuming the key is a single character
                }

                // Read the rest of the cells (values)
                while (std::getline(lineStream, cell, ',')) {
                    values.push_back(std::stoi(cell));
                }
                this->featureInfo[key] = values;
            }

            featureInfo_file.close();
        }
        
        void read_star_neighbors() {
            // read the star_neighbors from a csv file
            std::ifstream star_neighbors_file("required_files/border_starNeighbors/starNeighbors.csv");

            if (!star_neighbors_file.is_open()) {
                std::cerr << "Error opening file" << std::endl;
            }

            std::string line1;
            std::getline(star_neighbors_file, line1);  // Skip the first row (header)

            while (std::getline(star_neighbors_file, line1)) {
                std::stringstream ss(line1);
                int key;
                // Skip leading whitespace characters
                ss >> std::ws;
                ss >> key;
                ss.ignore(); // Ignore the comma

                // Read the list of ints
                int value;
                std::vector<int> values;
                while (ss >> value) {
                    values.push_back(value);
                    ss.ignore(); // Ignore the space
                }

                // Store the data in the map
                star_neighbors[key] = values;
            }

            star_neighbors_file.close();
        }
        
        std::vector<std::vector<std::string>> generate_size2_combos() {
            std::vector<std::string> features;
            for (const auto &entry : this->featureInfo) {
                features.push_back(entry.first);
            }

            std::vector<std::vector<std::string>> size2_candidatePatterns;
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
        
        std::vector<std::vector<std::string>> degree2Processing(
            std::vector<std::vector<std::string>> candidatePatterns, int candidatePatterns_size, 
            double prevalence_threshold) {   
            std::vector<std::vector<std::string>> size2_patterns;
            
            for (const auto& coloc : candidatePatterns) {
                std::string first_feature = coloc[0];
                std::string second_feature = coloc[1];
                
                // determine first feature start and end
                auto feature1 = this->featureInfo.find(first_feature);
                std::vector<int> values1 = feature1->second;
                int first_feature_start = values1[1];
                int first_feature_end = values1[2];
                // determine second feature start and end
                auto feature2 = this->featureInfo.find(second_feature);
                std::vector<int> values2 = feature2->second;
                int second_feature_start = values2[1];
                int second_feature_end = values2[2];
                
                std::vector<std::string> coloc_key = {first_feature, second_feature};
                instance_table[coloc_key] = {};
                this->hashmap[coloc_key] = {};
                this->hashmap[coloc_key][first_feature] = {};
                this->hashmap[coloc_key][second_feature] = {};
                
                for (int index = first_feature_start; index <= first_feature_end; index++) {
                    auto star_neighbor_it = this->star_neighbors.find(index);
                    std::vector<int> neighbors = 
                        findNeighborsInRange(star_neighbor_it->second, second_feature_start,
                                            second_feature_end);
                    if (!neighbors.empty()) {
                        std::vector<int> index_tuple = {index};
                        this->instance_table[coloc_key][index_tuple] = neighbors;
                        this->hashmap[coloc_key][first_feature].insert(index);
                        for (int neighbor : neighbors) {
                            this->hashmap[coloc_key][second_feature].insert(neighbor);
                        }
                    }
                }
                
                double pr_first_feature = 
                    static_cast<double>(this->hashmap[coloc_key][first_feature].size()) /
                    this->featureInfo[first_feature][0];
                double pr_second_feature = 
                    static_cast<double>(this->hashmap[coloc_key][second_feature].size()) /
                    this->featureInfo[second_feature][0];
                double PI = 0.0;

                if (pr_first_feature < pr_second_feature) {
                    PI = pr_first_feature;
                } else {
                    PI = pr_second_feature;
                }

                if (PI >= prevalence_threshold) {
                       size2_patterns.push_back(coloc_key);
                }
            }
            
            std::cout << "Degree 2 Prevalent Patterns for Border " + 
                std::to_string(this->border_id) + ":" << std::endl;
            for (auto i : size2_patterns) {
                std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
            }
            return size2_patterns;
        }
    };
    std::vector<Border> borders;

    // this class holds all information pertaining to the entire region
    class Region {
    public:
        std::map<std::vector<std::string>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<std::string>, std::map<std::string, std::set<int>>> hashmap;
        
        std::vector<std::vector<std::string>> generate_size2_combos() {
            std::set<std::vector<std::string>> patterns;
            for (const auto& entry : this->hashmap) {
                patterns.insert(entry.first);
            }
            std::vector<std::vector<std::string>> size2_candidatePatterns(patterns.begin(), patterns.end());
            return size2_candidatePatterns;
        }
        
        std::vector<std::vector<std::string>> degree2Processing(std::vector<std::vector<std::string>> candidatePatterns, 
                                                                int candidatePatterns_size, double prevalence_threshold,
                                                                int number_subregions) {
            std::vector<std::vector<std::string>> prevalent;
            
            for (const auto& coloc : candidatePatterns) {
                std::string first_feature = coloc[0];
                std::string second_feature = coloc[1];
                int total_first_feature = 0;
                int total_second_feature = 0;
                
                for (int i = 0; i < number_subregions; i++) {
                    auto feature1 = subregions[i].featureInfo.find(first_feature);
                    std::vector<int> values1 = feature1->second;
                    total_first_feature += values1[0];
                    
                    auto feature2 = subregions[i].featureInfo.find(second_feature);
                    std::vector<int> values2 = feature2->second;
                    total_second_feature += values2[0];
                }
                
                double pr_first_feature = static_cast<double>(this->hashmap[coloc][first_feature].size()) / total_first_feature;
                double pr_second_feature = static_cast<double>(this->hashmap[coloc][second_feature].size()) / total_second_feature;
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
            std::ofstream size2_file("size2_patterns_entire_region.txt");
            for (auto i : prevalent) {
                std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
                size2_file << "(" << i[0] << ", " << i[1] << ")" << std::endl;
            }
            size2_file.close();
            return prevalent;
        }
        
        void generateCombinations(const std::vector<std::string>& features, int degree, 
                                  std::vector<std::vector<std::string>>& result, 
                                  std::vector<std::string>& current, int start) {
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
        bool allSubpatternsInPrevalent(const std::vector<std::string>& pattern, 
                                       const std::set<std::vector<std::string>>& prevalentPatterns, 
                                       int degree) {
            std::vector<std::vector<std::string>> subpatterns;
            std::vector<std::string> current;
            generateCombinations(pattern, degree - 1, subpatterns, current, 0);

            for (const auto& subpattern : subpatterns) {
                if (prevalentPatterns.find(subpattern) == prevalentPatterns.end()) {
                    return false;
                }
            }
            return true;
        }

        std::vector<std::vector<std::string>> getCandidatePatterns(
            const std::vector<std::vector<std::string>>& prevalentPattern, int degree, std::vector<std::string> features) {
            std::set<std::vector<std::string>> prevalentPatterns(prevalentPattern.begin(), prevalentPattern.end());
            
            std::vector<std::vector<std::string>> _candidatePatterns;
            std::vector<std::vector<std::string>> _patterns;
            std::vector<std::string> current;

            generateCombinations(features, degree, _patterns, current, 0);

            for (const auto& pattern : _patterns) {
                if (allSubpatternsInPrevalent(pattern, prevalentPatterns, degree)) {
                    _candidatePatterns.push_back(pattern);
                }
            }

            return _candidatePatterns;
        }
        
        std::vector<int> findNeighborsInRange(const std::vector<int>& arr, int x, int y) {
            auto start_it = std::lower_bound(arr.begin(), arr.end(), x);
            auto end_it = std::upper_bound(arr.begin(), arr.end(), y);
            return std::vector<int>(start_it, end_it);
        }
        
        std::vector<std::vector<std::string>> colocationGeneral(std::vector<std::vector<std::string>> 
                                                                candidatePatterns, 
                                                                int candidatePatterns_size, 
                                                                double prevalence_threshold,
                                                                int degree, int number_subregions) {
            std::vector<std::vector<std::string>> prevalent;
            
            for (const auto& currentPattern : candidatePatterns) {
                std::vector<std::string> basePattern;
                for (int j = 0; j < degree - 1; j++) {
                    basePattern.push_back(currentPattern[j]);
                }
                std::string lastFeature = currentPattern[degree - 1];
                // Add entries to instance_table and hashmap
                this->instance_table[currentPattern] = {};
                this->hashmap[currentPattern] = {};
                
                // Initialize sets for each element in currentPattern in hashmap
                for (const auto& f : currentPattern) {
                    this->hashmap[currentPattern][f] = {};
                }
                
                auto colocTableIt = this->instance_table.find(basePattern);
                std::map<std::vector<int>, std::vector<int>>& colocTable = colocTableIt->second;
                
                for (const auto& entry : colocTable) {
                    const std::vector<int>& key = entry.first;
                    std::set<int> commonLastNeighbors;
                    for (int instanceID : key) {
                        std::set<int> lastFeatureNeighbors;
                        for (auto subregion : subregions) {
                            int lastFeatureStart = subregion.featureInfo[lastFeature][1];
                            int lastFeatureEnd = subregion.featureInfo[lastFeature][2];
                            
                            auto subregion_star_neighbor_it = subregion.star_neighbors.find(instanceID);
                            if (subregion_star_neighbor_it != subregion.star_neighbors.end()) {
                                std::vector<int> subregion_star_neighbor_list(subregion_star_neighbor_it->second.begin(), 
                                                   subregion_star_neighbor_it->second.end());
                                std::vector<int> neighborsInRange_vector = this->findNeighborsInRange(subregion_star_neighbor_list, lastFeatureStart, 
                                                                                      lastFeatureEnd);
                                std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());
                                lastFeatureNeighbors.insert(neighborsInRange.begin(), neighborsInRange.end());
                            }
                            auto border_star_neighbor_it = borders[0].star_neighbors.find(instanceID);
                            if (border_star_neighbor_it != borders[0].star_neighbors.end()) {
                                std::vector<int> border_star_neighbor_list(border_star_neighbor_it->second.begin(), 
                                                   border_star_neighbor_it->second.end());
                                std::vector<int> neighborsInRange_vector = this->findNeighborsInRange(border_star_neighbor_list, lastFeatureStart, 
                                                                                      lastFeatureEnd);
                                std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());
                                lastFeatureNeighbors.insert(neighborsInRange.begin(), neighborsInRange.end());
                            }
                        }
                        if (commonLastNeighbors.empty()) {
                            commonLastNeighbors = lastFeatureNeighbors;
                        } else {
                            std::set<int> intersection;
                            std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                                  lastFeatureNeighbors.begin(), lastFeatureNeighbors.end(),
                                                  std::inserter(intersection, intersection.begin()));
                            commonLastNeighbors = intersection;
                        }
                    }
                    for (int n : colocTable[key]) {
                        std::set<int> all_neighbors;
                        for (auto subregion : subregions) {
                            int lastFeatureStart = subregion.featureInfo[lastFeature][1];
                            int lastFeatureEnd = subregion.featureInfo[lastFeature][2];
                            std::set<int> neighbors_subregion;
                            auto subregion_star_neighbor_it = subregion.star_neighbors.find(n);
                            if (subregion_star_neighbor_it != subregion.star_neighbors.end()) {
                                std::vector<int> subregion_star_neighbor_list(subregion_star_neighbor_it->second.begin(), 
                                                   subregion_star_neighbor_it->second.end());
                                std::vector<int> neighborsInRange_vector = this->findNeighborsInRange(subregion_star_neighbor_list, lastFeatureStart, 
                                                                                      lastFeatureEnd);
                                std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());
                                std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                                     neighborsInRange.begin(), neighborsInRange.end(),
                                                     std::inserter(neighbors_subregion, neighbors_subregion.begin()));
                            }
                            std::set<int> neighbors_border;
                            auto border_star_neighbor_it = borders[0].star_neighbors.find(n);
                            if (border_star_neighbor_it != borders[0].star_neighbors.end()) {
                                std::vector<int> border_star_neighbor_list(border_star_neighbor_it->second.begin(), 
                                                   border_star_neighbor_it->second.end());
                                std::vector<int> neighborsInRange_vector = this->findNeighborsInRange(border_star_neighbor_list, lastFeatureStart, 
                                                                                      lastFeatureEnd);
                                std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());
                                std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                                     neighborsInRange.begin(), neighborsInRange.end(),
                                                     std::inserter(neighbors_border, neighbors_border.begin()));
                            }
                            all_neighbors.insert(neighbors_subregion.begin(), neighbors_subregion.end());
                            all_neighbors.insert(neighbors_border.begin(), neighbors_border.end());
                        }
                        if (!all_neighbors.empty()) {
                            std::vector<int> new_key = key;
                            new_key.push_back(n);
                            std::vector<int> intersectionVec(all_neighbors.begin(), all_neighbors.end());
                            this->instance_table[currentPattern][new_key] = intersectionVec;

                            for (size_t k = 0; k < new_key.size(); ++k) {
                                this->hashmap[currentPattern][currentPattern[k]].insert(new_key[k]);
                            }
                            this->hashmap[currentPattern][lastFeature].insert(all_neighbors.begin(), all_neighbors.end());
                        }
                    }
                }
                std::vector<double> pr;
                    for (int m = 0; m < degree; ++m) {
                        std::string f = currentPattern[m];
                        int total_count = 0;
                        for (auto subregion : subregions) {
                            total_count += subregion.featureInfo[f][0];
                        }
                        double ratio = static_cast<double>(this->hashmap[currentPattern][f].size()) 
                            / total_count;
                        pr.push_back(ratio);
                    }
                    double PI = *std::min_element(pr.begin(), pr.end());
                    if (PI >= prevalence_threshold) {
                        prevalent.push_back(currentPattern);
                    }
            }
            std::cout << "Degree " << degree << " Prevalent Patterns for Entire Region:" +
                std::to_string(number_subregions) + ":"<< std::endl;
            std::ofstream patterns_file("size" + std::to_string(degree) + "_patterns_entire_region.txt");
            for (auto& patternVec : prevalent) {
                std::cout << "(";
                patterns_file << "(";
                for (size_t i = 0; i < patternVec.size(); i++) {
                    std::cout << patternVec[i];
                    patterns_file << patternVec[i];
                    if (i < patternVec.size() - 1) {
                        std::cout << ", ";
                        patterns_file << ", ";
                    }
                }
                std::cout << ")" << std::endl;
                patterns_file << ")" << std::endl;
            }
            patterns_file.close();
            return prevalent;
        }
    };
    Region region;

    // processing for all the sub-regions
    void subregion_main(int number_subregions, double prevalence_threshold) {
        for (int i = 0; i < number_subregions; i++) {
            subregions.push_back(SubRegion(i));
        }
        
        for (int i = 0; i < number_subregions; i++) {
            subregions[i].read_featureInfo();
            subregions[i].read_star_neighbors();
            std::vector<std::vector<std::string>> size2_candidatePatterns = 
                subregions[i].generate_size2_combos();
            subregions[i].size2_patterns = 
                subregions[i].degree2Processing(size2_candidatePatterns,
                                                size2_candidatePatterns.size(), 
                                                prevalence_threshold, i);
            
            int degree = 3;
            std::vector<std::vector<std::string>> candidatePatterns =
                subregions[i].getCandidatePatterns(subregions[i].size2_patterns, degree);
            while (!candidatePatterns.empty()) {
                std::vector<std::vector<std::string>> prevalent_patterns = 
                subregions[i].colocationGeneral(candidatePatterns, candidatePatterns.size(), 
                           prevalence_threshold, degree, i);
                degree += 1;
                if (prevalent_patterns.size() == 0) {
                    break;
                }
                candidatePatterns = subregions[i].getCandidatePatterns(prevalent_patterns,
                                                                       degree);
            }
        }
    }

    // processing for the border region
    void border_main(int number_borders, double prevalence_threshold) {
        for (int i = 0; i < number_borders; i++) {
            borders.push_back(Border(i));
        }
        
        for (int i = 0; i < number_borders; i ++) {
            borders[i].read_featureInfo();
            borders[i].read_star_neighbors();
            std::vector<std::vector<std::string>> size2_candidatePatterns = 
                borders[i].generate_size2_combos();
            borders[i].size2_patterns = 
                borders[i].degree2Processing(size2_candidatePatterns, 
                                             size2_candidatePatterns.size(), 
                                             prevalence_threshold);
        } 
    } 

    void update_border_info(int* ids, int ids_size, int i) {
        // update the border hashmap so that it has the original IDS not the indicies
        std::map<std::vector<std::string>, std::map<std::string, std::set<int>>> temp_hash;
        for (const auto& outer_pair : borders[i].hashmap) {
            const std::vector<std::string>& outer_key = outer_pair.first;
            const auto& inner_map = outer_pair.second;
            std::map<std::string, std::set<int>> temp_inner_map;
            for (const auto& inner_pair : inner_map) {
                const std::string& inner_key = inner_pair.first;
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
        std::map<std::vector<std::string>, std::map<std::vector<int>, std::vector<int>>> temp_instance_table;

        for (const auto& outer_pair : borders[i].instance_table) {
            const std::vector<std::string>& outer_key = outer_pair.first;
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
        
        for (const auto& entry : borders[i].star_neighbors) {
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
        borders[i].star_neighbors = temp_star_neighbors;        
    }

    // combine the hashmaps of the subregions and border region
    void combine_hashmaps(int number_subregions, int number_borders) {
        std::set<std::vector<std::string>> keys_needed;
        
        // determine which keys from the subregions and borders need to be combined to form the regional hashmap
        for (int i = 0; i < number_subregions; i++) {
            std::set<std::vector<std::string>> temp_vect(subregions[i].size2_patterns.begin(), 
                                                         subregions[i].size2_patterns.end());  
            if (i == 0) {
                keys_needed = temp_vect;
            } else {
                std::set<std::vector<std::string>> temp_result;
                std::set_union(keys_needed.begin(), keys_needed.end(),
                               temp_vect.begin(), temp_vect.end(),
                               std::inserter(temp_result, temp_result.begin()));
                keys_needed = temp_result;
            }
        }
        
        for (int i = 0; i < number_borders; i++) {
            std::set<std::vector<std::string>> temp_vect(borders[i].size2_patterns.begin(), 
                                                         borders[i].size2_patterns.end());
            std::set<std::vector<std::string>> temp_result;
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

    // combine the instance tables of the subregions and the border region
    void combine_instance_tables(int number_subregions, int number_borders) {
        std::set<std::vector<std::string>> keys_needed;
        
        // determine which keys from the subregions and borders need to be combined to form the regional hashmap
        for (int i = 0; i < number_subregions; i++) {
            std::set<std::vector<std::string>> temp_vect(subregions[i].size2_patterns.begin(), 
                                                         subregions[i].size2_patterns.end());  
            if (i == 0) {
                keys_needed = temp_vect;
            } else {
                std::set<std::vector<std::string>> temp_result;
                std::set_union(keys_needed.begin(), keys_needed.end(),
                               temp_vect.begin(), temp_vect.end(),
                               std::inserter(temp_result, temp_result.begin()));
                keys_needed = temp_result;
            }
        }
        
        for (int i = 0; i < number_borders; i++) {
            std::set<std::vector<std::string>> temp_vect(borders[i].size2_patterns.begin(), 
                                                         borders[i].size2_patterns.end());
            std::set<std::vector<std::string>> temp_result;
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

    // processing for the entire region
    void region_main(int number_subregions, double prevalence_threshold, char** features_ptr, int features_size) {
        std::vector<std::string> features(features_ptr, features_ptr + features_size);
        std::vector<std::vector<std::string>> size2_candidatePatterns = region.generate_size2_combos();
        std::vector<std::vector<std::string>> prevalentPatterns = region.degree2Processing(size2_candidatePatterns, 
                                                                size2_candidatePatterns.size(), prevalence_threshold,
                                                                                           number_subregions);
        int degree = 3;
        std::vector<std::vector<std::string>> candidatePatterns = region.getCandidatePatterns(prevalentPatterns,
                                                                                              degree, features);
        while (!candidatePatterns.empty()) {
            std::vector<std::vector<std::string>> prevalent_patterns = region.colocationGeneral(candidatePatterns, 
                                                                    candidatePatterns.size(), prevalence_threshold, 
                                                                    degree, number_subregions);
            degree += 1;
            if (prevalent_patterns.size() == 0) {
                break;
            }
            candidatePatterns = region.getCandidatePatterns(prevalentPatterns, degree, features);
        }
    }
}
