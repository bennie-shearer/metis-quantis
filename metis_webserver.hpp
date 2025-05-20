#ifndef METIS_WEBSERVER_HPP
#define METIS_WEBSERVER_HPP

#include <string>
#include <functional>
#include <map>
#include <thread>

// Windows specific socket types
#ifdef _WIN32
    #include <winsock2.h>
    typedef SOCKET socket_t;
#else
typedef int socket_t;
#endif

namespace simple_web {
    class Request {
    public:
        std::string path;
        std::string method;
        std::string body;
        // Add other properties as needed
    };

    class Response {
    public:
        void send(const std::string& content) {
            // Implementation
        }

        void setStatus(int status) {
            // Implementation
        }

        void setHeader(const std::string& name, const std::string& value) {
            // Implementation
        }
        // Add other methods as needed
    };

    class SimpleWebServer {
    public:
        SimpleWebServer(int port);
        ~SimpleWebServer();

        // Method to register API endpoints
        void get(const std::string& path, std::function<void(const Request&, Response&)> handler);
        void post(const std::string& path, std::function<void(const Request&, Response&)> handler);

        // Method to set the static files path
        void setStaticFilesPath(const std::string& path);

        // Server control methods
        bool start();
        void stop();
        void run();
        bool isRunning() const;

    private:
        int port_;
        bool running_;
        socket_t serverSocket_;
        std::thread serverThread_;
        std::string staticFilesPath_;

        void serverLoop();
        std::string readFile(const std::string& filePath);
        std::string getContentType(const std::string& filePath);
    };
}

#endif // METIS_WEBSERVER_HPP