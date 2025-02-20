#ifndef LOGGER_H
#define LOGGER_H

#include <sstream>
#include <fstream>
#include <chrono>
#include <iostream>
class Logger {
public:
    static void init(const std::string& filename) {
        logFile.open(filename);
    }
    
    template<typename T>
    static void log(const T& message, bool console = true) {
        auto now = std::chrono::system_clock::now();
        std::stringstream ss;
        ss << "[" << std::chrono::system_clock::to_time_t(now) << "] " << message;
        
        if(console) std::cout << ss.str() << std::endl;
        if(logFile.is_open()) logFile << ss.str() << std::endl;
    }

private:
    static std::ofstream logFile;
};

#endif
