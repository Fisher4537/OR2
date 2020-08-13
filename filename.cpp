#include <string>
#include <iostream>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

extern "C" {

    char** get_all_data_filenames(int folders) {

        /* folders:
            0 = light
            1 = average
            2 = heavy
            3 = light + average
            4 = light + average + heavy
        */


        folders = folders == 0 ? 29 : folders == 1 ? 16 : folders == 2 ? 43 : folders == 3 ? 45 : 88;

        char** files = (char**)calloc(folders, sizeof(char*));
        for(int i = 0; i < folders; i++)
            files[i] = (char*)calloc(25, sizeof(char));

        std::string path = "\data_light";
        int idx = 0;
        for (const auto& entry : fs::directory_iterator(path)) {
            std::cout << entry.path() << std::endl;
            idx++;
        }
    }
}