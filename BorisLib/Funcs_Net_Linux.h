#pragma once

#include "BorisLib_Config.h"
#if OPERATING_SYSTEM == OS_LIN

#include <netdb.h>
#include <arpa/inet.h>
#include <string.h>

#include "Funcs_Conv.h"
#include "LinSocks.h"

//----------------------------------------------------------------------------------------------------------------------

//open a web page with given url in default system browser
inline bool open_web_page(std::string url)
{
	if (url.find(std::string("http://")) == std::string::npos && url.find(std::string("https://")) == std::string::npos)
		url = std::string("http://") + url;

	std::string argument = std::string("xdg-open ") + url;
	int retval = system(argument.c_str());
	
	//if you don't use the return value of system call g++ will spit out an annoying warning (which could be turned off, but better to get retval)
	return retval >= 0;
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

		struct addrinfo *result = nullptr;
		struct addrinfo hints;

		int errcode;
		char addrstr[100];
		void* ptr;
	
		memset(&hints, 0, sizeof(hints));
		hints.ai_family = PF_UNSPEC;
  		hints.ai_socktype = SOCK_STREAM;
  		hints.ai_flags |= AI_CANONNAME;

  		errcode = getaddrinfo(domain_name.c_str(), NULL, &hints, &result);
  		if (errcode != 0) {
		  
		  //failed
		  freeaddrinfo(result);
		  return false;
	  	}

  		while (result) {
     		
			switch (result->ai_family) {
     			
				case AF_INET:
				{
        			ptr = &((struct sockaddr_in*) result->ai_addr)->sin_addr;
        			inet_ntop(result->ai_family, ptr, addrstr, 100);

					ip_address = std::string(addrstr);
				}
        		break;

     			case AF_INET6:
				{
        			ptr = &((struct sockaddr_in6*) result->ai_addr)->sin6_addr;
        			inet_ntop(result->ai_family, ptr, addrstr, 100);

					ip_address = std::string(addrstr);
				}
        		break;          
     		}

     		result = result->ai_next;
  		}

		freeaddrinfo(result);
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