#include "metis_webserver.hpp"
#include <iostream>
#include <thread>
#include <chrono>
#include <string>
#include <map>
#include <functional>
#include <cstring>
#include <fstream>
#include <sstream>
#include <filesystem>

// Include socket libraries for Windows
#ifdef _WIN32
    #include <winsock2.h>
    #include <ws2tcpip.h>
    #pragma comment(lib, "ws2_32.lib")
#else
    #include <sys/socket.h>
    #include <netinet/in.h>
    #include <unistd.h>
    #include <arpa/inet.h>
    #include <cerrno>
#endif

namespace simple_web {

    SimpleWebServer::SimpleWebServer(int port)
        : port_(port), running_(false), staticFilesPath_("html") {

        #ifdef _WIN32
            // Initialize Winsock on Windows
            WSADATA wsaData;
            if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) {
                std::cerr << "Failed to initialize Winsock" << std::endl;
                return;
            }
        #endif

        // Initialize any other members
    }

    SimpleWebServer::~SimpleWebServer() {
        if (running_) {
            stop();
        }

        #ifdef _WIN32
            WSACleanup();
        #endif
    }

    bool SimpleWebServer::start() {
        #ifdef _WIN32
            // Create a socket
            serverSocket_ = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
            if (serverSocket_ == INVALID_SOCKET) {
                std::cerr << "Error opening socket: " << WSAGetLastError() << std::endl;
                return false;
            }

            // Set socket options
            BOOL opt = TRUE;
            if (setsockopt(serverSocket_, SOL_SOCKET, SO_REUSEADDR, (char*)&opt, sizeof(opt)) == SOCKET_ERROR) {
                std::cerr << "Error setting socket options: " << WSAGetLastError() << std::endl;
                closesocket(serverSocket_);
                return false;
            }

            // Bind the socket
            struct sockaddr_in serverAddr;
            serverAddr.sin_family = AF_INET;
            serverAddr.sin_addr.s_addr = INADDR_ANY;
            serverAddr.sin_port = htons(port_);

            if (bind(serverSocket_, (struct sockaddr*)&serverAddr, sizeof(serverAddr)) == SOCKET_ERROR) {
                std::cerr << "Error on binding: " << WSAGetLastError() << std::endl;
                closesocket(serverSocket_);
                return false;
            }

            // Start listening for connections
            if (listen(serverSocket_, SOMAXCONN) == SOCKET_ERROR) {
                std::cerr << "Error on listen: " << WSAGetLastError() << std::endl;
                closesocket(serverSocket_);
                return false;
            }

            // Start the server thread
            running_ = true;
            serverThread_ = std::thread(&SimpleWebServer::serverLoop, this);

            std::cout << "Server started on all interfaces on port " << port_ << std::endl;
            return true;
        #else
            // Unix/Linux implementation
            serverSocket_ = socket(AF_INET, SOCK_STREAM, 0);
            if (serverSocket_ < 0) {
                std::cerr << "Error opening socket" << std::endl;
                return false;
            }

            // Set socket options
            int opt = 1;
            if (setsockopt(serverSocket_, SOL_SOCKET, SO_REUSEADDR, (const char*)&opt, sizeof(opt)) < 0) {
                std::cerr << "Error setting socket options" << std::endl;
                close(serverSocket_);
                return false;
            }

            // Bind the socket
            struct sockaddr_in serverAddr;
            memset(&serverAddr, 0, sizeof(serverAddr));
            serverAddr.sin_family = AF_INET;
            serverAddr.sin_addr.s_addr = INADDR_ANY;
            serverAddr.sin_port = htons(port_);

            if (bind(serverSocket_, (struct sockaddr*)&serverAddr, sizeof(serverAddr)) < 0) {
                std::cerr << "Error on binding: " << strerror(errno) << std::endl;
                close(serverSocket_);
                return false;
            }

            // Start listening for connections
            if (listen(serverSocket_, 5) < 0) {
                std::cerr << "Error on listen: " << strerror(errno) << std::endl;
                close(serverSocket_);
                return false;
            }

            // Start the server thread
            running_ = true;
            serverThread_ = std::thread(&SimpleWebServer::serverLoop, this);

            std::cout << "Server started on all interfaces on port " << port_ << std::endl;
            return true;
        #endif
    }

    void SimpleWebServer::stop() {
        running_ = false;

        if (serverThread_.joinable()) {
            serverThread_.join();
        }

        #ifdef _WIN32
            closesocket(serverSocket_);
        #else
            close(serverSocket_);
        #endif

        std::cout << "Server stopped" << std::endl;
    }

    void SimpleWebServer::run() {
        start();

        // Keep the main thread alive
        while (running_) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
    }

    void SimpleWebServer::get(const std::string& path, std::function<void(const Request&, Response&)> handler) {
        // Register a GET handler
        std::cout << "Registered GET handler for path: " << path << std::endl;
        // Implementation would store the handler for the path
    }

    void SimpleWebServer::post(const std::string& path, std::function<void(const Request&, Response&)> handler) {
        // Register a POST handler
        std::cout << "Registered POST handler for path: " << path << std::endl;
        // Implementation would store the handler for the path
    }

    void SimpleWebServer::setStaticFilesPath(const std::string& path) {
        staticFilesPath_ = path;
    }

    std::string SimpleWebServer::readFile(const std::string& filePath) {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filePath << std::endl;
            return "";
        }

        std::stringstream buffer;
        buffer << file.rdbuf();
        return buffer.str();
    }

    std::string SimpleWebServer::getContentType(const std::string& filePath) {
        // Get the file extension
        std::string extension = "";
        size_t dotPos = filePath.rfind('.');
        if (dotPos != std::string::npos) {
            extension = filePath.substr(dotPos + 1);
        }

        // Convert to lowercase
        for (auto& c : extension) {
            c = std::tolower(c);
        }

        // Return the appropriate content type based on the extension
        if (extension == "html" || extension == "htm") {
            return "text/html";
        } else if (extension == "css") {
            return "text/css";
        } else if (extension == "js") {
            return "application/javascript";
        } else if (extension == "json") {
            return "application/json";
        } else if (extension == "png") {
            return "image/png";
        } else if (extension == "jpg" || extension == "jpeg") {
            return "image/jpeg";
        } else if (extension == "gif") {
            return "image/gif";
        } else if (extension == "svg") {
            return "image/svg+xml";
        } else if (extension == "ico") {
            return "image/x-icon";
        } else if (extension == "txt") {
            return "text/plain";
        } else if (extension == "pdf") {
            return "application/pdf";
        } else {
            return "application/octet-stream";
        }
    }

    void SimpleWebServer::serverLoop() {
        std::cout << "Server loop started, waiting for connections..." << std::endl;

        while (running_) {
            #ifdef _WIN32
                // Accept connections on Windows
                struct sockaddr_in clientAddr;
                int clientLen = sizeof(clientAddr);

                std::cout << "Waiting to accept a connection..." << std::endl;
                SOCKET clientSocket = accept(serverSocket_, (struct sockaddr*)&clientAddr, &clientLen);

                if (clientSocket == INVALID_SOCKET) {
                    std::cerr << "Error on accept: " << WSAGetLastError() << std::endl;
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    continue;
                }

                // Handle the connection
                std::cout << "Connection accepted from "
                          << inet_ntoa(clientAddr.sin_addr) << ":"
                          << ntohs(clientAddr.sin_port) << std::endl;

                // Buffer to receive client data
                char buffer[4096] = {0};
                int bytesReceived = recv(clientSocket, buffer, sizeof(buffer) - 1, 0);

                if (bytesReceived > 0) {
                    // Extract the requested path from the HTTP request
                    std::string request(buffer);
                    size_t pathStart = request.find(' ') + 1;
                    size_t pathEnd = request.find(' ', pathStart);
                    std::string path = request.substr(pathStart, pathEnd - pathStart);

                    std::cout << "Requested path: " << path << std::endl;

                    // Check if path is root or specific file
                    if (path == "/" || path == "/index.html") {
                        path = "/index.html"; // Default to index.html
                    }

                    // Serve the file if it exists
                    std::string filePath = staticFilesPath_ + path;

                    // Normalize the path to avoid directory traversal attacks
                    std::filesystem::path normalizedPath = std::filesystem::weakly_canonical(std::filesystem::path(filePath));
                    std::string normalizedPathStr = normalizedPath.string();

                    // Check if the file exists
                    if (std::filesystem::exists(normalizedPath) && !std::filesystem::is_directory(normalizedPath)) {
                        // Read the file content
                        std::string fileContent = readFile(normalizedPathStr);

                        if (!fileContent.empty()) {
                            // Determine the content type
                            std::string contentType = getContentType(normalizedPathStr);

                            // Prepare the HTTP response
                            std::string httpResponse =
                                "HTTP/1.1 200 OK\r\n"
                                "Content-Type: " + contentType + "; charset=UTF-8\r\n"
                                "Content-Length: " + std::to_string(fileContent.length()) + "\r\n"
                                "Connection: close\r\n"
                                "\r\n";

                            // Send headers first
                            send(clientSocket, httpResponse.c_str(), httpResponse.length(), 0);

                            // Then send the file content
                            send(clientSocket, fileContent.c_str(), fileContent.length(), 0);

                            std::cout << "Served file: " << normalizedPathStr << std::endl;
                        } else {
                            // File could not be read
                            std::string errorResponse =
                                "HTTP/1.1 500 Internal Server Error\r\n"
                                "Content-Type: text/html; charset=UTF-8\r\n"
                                "Connection: close\r\n"
                                "\r\n"
                                "<html><body><h1>500 Internal Server Error</h1><p>Could not read file: " + path + "</p></body></html>";

                            send(clientSocket, errorResponse.c_str(), errorResponse.length(), 0);
                        }
                    } else {
                        // File not found, send 404 response
                        std::string notFoundResponse =
                            "HTTP/1.1 404 Not Found\r\n"
                            "Content-Type: text/html; charset=UTF-8\r\n"
                            "Connection: close\r\n"
                            "\r\n"
                            "<html><body><h1>404 Not Found</h1><p>The requested file was not found: " + path + "</p></body></html>";

                        send(clientSocket, notFoundResponse.c_str(), notFoundResponse.length(), 0);
                    }
                }

                // Close the client socket
                closesocket(clientSocket);
            #else
                // Accept connections on Unix/Linux
                struct sockaddr_in clientAddr;
                socklen_t clientLen = sizeof(clientAddr);

                std::cout << "Waiting to accept a connection..." << std::endl;
                int clientSocket = accept(serverSocket_, (struct sockaddr*)&clientAddr, &clientLen);

                if (clientSocket < 0) {
                    std::cerr << "Error on accept: " << strerror(errno) << std::endl;
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    continue;
                }

                // Handle the connection
                std::cout << "Connection accepted from "
                          << inet_ntoa(clientAddr.sin_addr) << ":"
                          << ntohs(clientAddr.sin_port) << std::endl;

                // Buffer to receive client data
                char buffer[4096] = {0};
                int bytesReceived = recv(clientSocket, buffer, sizeof(buffer) - 1, 0);

                if (bytesReceived > 0) {
                    // Extract the requested path from the HTTP request
                    std::string request(buffer);
                    size_t pathStart = request.find(' ') + 1;
                    size_t pathEnd = request.find(' ', pathStart);
                    std::string path = request.substr(pathStart, pathEnd - pathStart);

                    std::cout << "Requested path: " << path << std::endl;

                    // Check if path is root or specific file
                    if (path == "/" || path == "/index.html") {
                        path = "/index.html"; // Default to index.html
                    }

                    // Serve the file if it exists
                    std::string filePath = staticFilesPath_ + path;

                    // Normalize the path to avoid directory traversal attacks
                    std::filesystem::path normalizedPath = std::filesystem::weakly_canonical(std::filesystem::path(filePath));
                    std::string normalizedPathStr = normalizedPath.string();

                    // Check if the file exists
                    if (std::filesystem::exists(normalizedPath) && !std::filesystem::is_directory(normalizedPath)) {
                        // Read the file content
                        std::string fileContent = readFile(normalizedPathStr);

                        if (!fileContent.empty()) {
                            // Determine the content type
                            std::string contentType = getContentType(normalizedPathStr);

                            // Prepare the HTTP response
                            std::string httpResponse =
                                "HTTP/1.1 200 OK\r\n"
                                "Content-Type: " + contentType + "; charset=UTF-8\r\n"
                                "Content-Length: " + std::to_string(fileContent.length()) + "\r\n"
                                "Connection: close\r\n"
                                "\r\n";

                            // Send headers first
                            send(clientSocket, httpResponse.c_str(), httpResponse.length(), 0);

                            // Then send the file content
                            send(clientSocket, fileContent.c_str(), fileContent.length(), 0);

                            std::cout << "Served file: " << normalizedPathStr << std::endl;
                        } else {
                            // File could not be read
                            std::string errorResponse =
                                "HTTP/1.1 500 Internal Server Error\r\n"
                                "Content-Type: text/html; charset=UTF-8\r\n"
                                "Connection: close\r\n"
                                "\r\n"
                                "<html><body><h1>500 Internal Server Error</h1><p>Could not read file: " + path + "</p></body></html>";

                            send(clientSocket, errorResponse.c_str(), errorResponse.length(), 0);
                        }
                    } else {
                        // File not found, send 404 response
                        std::string notFoundResponse =
                            "HTTP/1.1 404 Not Found\r\n"
                            "Content-Type: text/html; charset=UTF-8\r\n"
                            "Connection: close\r\n"
                            "\r\n"
                            "<html><body><h1>404 Not Found</h1><p>The requested file was not found: " + path + "</p></body></html>";

                        send(clientSocket, notFoundResponse.c_str(), notFoundResponse.length(), 0);
                    }
                }

                // Close the client socket
                close(clientSocket);
            #endif

            // Slight delay to avoid excessive CPU usage
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }

    bool SimpleWebServer::isRunning() const {
        return running_;
    }
}