#include <iostream>
#include <vector>
#include <string>
#include <cstdio>

#ifdef _WIN32
    #include <process.h>
#else
    #include <unistd.h>
#endif

int main(int argc, char* argv[]) {
    // 1. Get the invocation path of the wrapper program itself
    std::string wrapper_path = argv[0];
    std::string dir_path = "";
    
    // 2. Find the last slash in the path (compatible with Windows '\' and Unix '/')
    size_t last_slash = wrapper_path.find_last_of("/\\");
    if (last_slash != std::string::npos) {
        // If found, extract the directory path, e.g., "/home/user/bin/" or "./"
        dir_path = wrapper_path.substr(0, last_slash + 1);
    }
    
    // 3. Look for the original USalign program in the extracted directory
#ifdef _WIN32
    std::string exe_name = dir_path + "USalign.exe";
#else
    std::string exe_name = dir_path + "USalign";
#endif

    std::vector<char*> args;

    // The first argument must be the executable name
    args.push_back(const_cast<char*>(exe_name.c_str()));

    // Pass all user-provided arguments exactly as they are
    for (int i = 1; i < argc; ++i) {
        args.push_back(argv[i]);
    }

    // Forcefully append the flexible alignment arguments at the end
    args.push_back(const_cast<char*>("-mm"));
    args.push_back(const_cast<char*>("7"));
    args.push_back(nullptr); // The argument array must be null-terminated

    // 4. Transfer execution control to the original USalign
#ifdef _WIN32
    intptr_t exit_code = _spawnvp(_P_WAIT, exe_name.c_str(), args.data());
    if (exit_code == -1) {
        std::cerr << "Execution failed! Please ensure " << exe_name 
                  << " is in the same directory as this program." << std::endl;
        return 1;
    }
    return static_cast<int>(exit_code);
#else
    // execvp replaces the current process with the target program
    execvp(exe_name.c_str(), args.data());
    
    // The following lines are only reached if execvp fails
    std::string error_msg = "Execution failed! Please ensure " + exe_name + " exists and is executable";
    perror(error_msg.c_str());
    return 1;
#endif
}