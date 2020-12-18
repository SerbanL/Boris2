#pragma once

#include "BorisLib_Config.h"
#if OPERATING_SYSTEM == OS_WIN

#include <shellapi.h>

#include "Funcs_Conv.h"
#include "Funcs_Conv_Windows.h"
#include "WinSocks.h"

//----------------------------------------------------------------------------------------------------------------------

//open a web page with given url in default system browser
inline bool open_web_page(std::string url)
{
	ShellExecute(0, 0, StringtoWCHARPointer(url), 0, 0, SW_SHOW);
	//see corresponding function in Funcs_Net_Linux.h for why a return is needed
	return true;
}

//----------------------------------------------------------------------------------------------------------------------

class HTTP
{
private:

	std::string ip_address;
	std::string domain_name;

	bool available = false;

private:

	//From given domain name and port (default 80) get the ip address as a std::string (set it in HTTP::ip_address); also save domain_name to HTTP::domain_name. Return false if failed.
	bool ip_address_from_hostname(std::string domain_name)
	{
		this->domain_name = domain_name;

		//from https://docs.microsoft.com/en-us/windows/desktop/api/ws2tcpip/nf-ws2tcpip-getaddrinfo

		//--------------------------------
		// Declare and initialize variables
		WSADATA wsaData;
		int iResult;
		INT iRetval;

		DWORD dwRetval;

		int i = 1;

		struct addrinfo *result = NULL;
		struct addrinfo *ptr = NULL;
		struct addrinfo hints;

		struct sockaddr_in  *sockaddr_ipv4;

		LPSOCKADDR sockaddr_ip;

		char ipstringbuffer[46];
		DWORD ipbufferlength = 46;

		// Initialize Winsock
		iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
		if (iResult != 0) {

			//failed
			return false;
		}

		//--------------------------------
		// Setup the hints address info structure
		// which is passed to the getaddrinfo() function
		ZeroMemory(&hints, sizeof(hints));
		hints.ai_family = AF_UNSPEC;
		hints.ai_socktype = SOCK_STREAM;
		hints.ai_protocol = IPPROTO_TCP;

		//--------------------------------
		// Call getaddrinfo(). If the call succeeds,
		// the result variable will hold a linked list
		// of addrinfo structures containing response
		// information
		dwRetval = getaddrinfo(domain_name.c_str(), "80", &hints, &result);
		if (dwRetval != 0) {

			//failed
			WSACleanup();
			return false;
		}

		// Retrieve each address and print out the hex BYTEs
		for (ptr = result; ptr != NULL; ptr = ptr->ai_next) {

			switch (ptr->ai_family) {

				//IPv4
			case AF_INET:
			{
				sockaddr_ipv4 = (struct sockaddr_in *) ptr->ai_addr;
				ip_address = std::string(inet_ntoa(sockaddr_ipv4->sin_addr));
			}
			break;

			//IPv6
			case AF_INET6:
			{
				// We use WSAAddressToString since it is supported on Windows XP and later
				sockaddr_ip = (LPSOCKADDR)ptr->ai_addr;
				// The buffer length is changed by each call to WSAAddresstoString
				// So we need to set it for each iteration through the loop for safety
				ipbufferlength = 46;
				iRetval = WSAAddressToString(sockaddr_ip, (DWORD)ptr->ai_addrlen, nullptr, StringtoWCHARPointer(std::string(ipstringbuffer)), &ipbufferlength);
				if (iRetval) {

					//error
					return false;
				}
				else {

					ip_address = std::string(ipstringbuffer);
				}
			}
			break;

			default:
				break;
			}
		}

		freeaddrinfo(result);
		WSACleanup();

		return true;
	}

public:

	HTTP(std::string domain_name)
	{
		available = ip_address_from_hostname(domain_name);
	}

	bool is_available(void) { return available; }

	std::string get_ip_address(void) { return ip_address; }

	//send an HTTP 1.1 POST request to a domainname/page (e.g. to www.boris-spintronics.uk/version.php, where domainname = www.boris-spintronics.uk and page = version.php) and obtain response
	//send POST parameters as a std::string using the usual POST syntax
	//return false if couldn't obtain response. If true then received message set in response
	bool http_post(std::string page, std::string parameters, std::string& response_message, int buffer_size = DEFAULT_BUFLEN)
	{
		//ip address must be set in constructor correctly
		if (!available) return false;

		NetSocksClient client(ip_address, std::string("80"), buffer_size);

		std::string message = std::string("POST /") + page + " HTTP/1.1\r\n";
		message += std::string("Host: ") + domain_name + "\r\n";
		message += std::string("Content-Type: application/x-www-form-urlencoded") + "\r\n";
		message += std::string("Content-Length: ");
		message += ToString(parameters.length()) + "\r\n";
		message += std::string("\r\n");
		message += parameters + "\r\n";

		std::string response = client.SendMessage_GetResponse(message);

		std::string OK_message = "HTTP/1.1 200 OK";

		if (response.length() > OK_message.length() && response.substr(0, 15) == OK_message) {

			//successful connection attempt and message received : response must start with "HTTP/1.1 200 OK"

			//message starts after "CRLFCRLF", i.e. "\r\n\r\n"
			size_t pos_start = response.find("\r\n\r\n");

			if (pos_start != std::string::npos) {

				//between pos_start + 4 and next "CRLF" we have the number of message characters in hex
				std::string message_no_header = response.substr(pos_start + 4);

				size_t pos_message_start = message_no_header.find("\r\n");

				if (pos_message_start != std::string::npos) {

					std::string message_length_hex = message_no_header.substr(0, pos_message_start);

					//convert from hex std::string to int
					std::stringstream ss;
					ss << std::hex << message_length_hex;
					int num_chars;
					ss >> num_chars;

					response_message = message_no_header.substr(pos_message_start + 2, num_chars);

					return true;
				}
			}
		}

		return false;
	}

	//send an HTTP 1.1 GET request to a domainname/page (e.g. to www.boris-spintronics.uk/version.php, where domainname = www.boris-spintronics.uk and page = version.php) and obtain response
	//return false if couldn't obtain response. If true then received message set in response
	bool http_get(std::string page, std::string& response_message)
	{
		//ip address must be set in constructor correctly
		if (!available) return false;

		NetSocksClient client(ip_address, std::string("80"));

		std::string message = std::string("GET /") + page + " HTTP/1.1\r\n";
		message += std::string("Host: ") + domain_name + "\r\n";
		message += std::string("Connection: Close\r\n\r\n");

		std::string response = client.SendMessage_GetResponse(message);

		std::string OK_message = "HTTP/1.1 200 OK";

		if (response.length() > OK_message.length() && response.substr(0, 15) == OK_message) {

			//successful connection attempt and message received : response must start with "HTTP/1.1 200 OK"

			//message starts after "CRLFCRLF", i.e. "\r\n\r\n"
			size_t pos_start = response.find("\r\n\r\n");

			if (pos_start != std::string::npos) {

				//between pos_start + 4 and next "CRLF" we have the number of message characters in hex
				std::string message_no_header = response.substr(pos_start + 4);

				size_t pos_message_start = message_no_header.find("\r\n");

				if (pos_message_start != std::string::npos) {

					std::string message_length_hex = message_no_header.substr(0, pos_message_start);

					//convert from hex std::string to int
					std::stringstream ss;
					ss << std::hex << message_length_hex;
					int num_chars;
					ss >> num_chars;

					response_message = message_no_header.substr(pos_message_start + 2, num_chars);

					return true;
				}
			}
		}

		return false;
	}
};

#endif