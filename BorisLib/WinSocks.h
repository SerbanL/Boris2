#pragma once

#include <Windows.h>
#include <Winsock2.h>
#include <ws2tcpip.h>
#include <string>
#include <vector>

#pragma comment (lib, "Ws2_32.lib")

#define DEFAULT_BUFLEN 1024
#define DEFAULT_PORT "1542"

#define RECVSLEEPMS		100				//Used for non-blocking calls to recv : after each call sleep this number of ms to avoid hogging too much CPU time 
										//Calling recv too often, or using blocking sockets i.e. blocking recv, will put a slight load on one of the CPU cores which can really degrade parallelized computations, irrespective of which thread what

//see https://msdn.microsoft.com/en-us/library/windows/desktop/ms737889(v=vs.85).aspx

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//SERVER

class WinSocks {

private:

    SOCKET ListenSocket;
    SOCKET ClientSocket;

	//receive buffer and length
    char *recvbuf;
    int recvbuflen;

	//flags
	bool clientConnected;
	bool listenSocketActive;

	//commands processed by Simulation may return parameters. These are placed in this std::vector to be sent to Client in winsocks thread
	std::vector<std::string> dataParams;

	//the last received message
	std::string last_message;

private:

	//Make non-blocking listen socket on DEFAULT_PORT to accept incoming connections - done in the constructor
	void MakeListenSocket(void);

public:

	WinSocks(int buffer_size = DEFAULT_BUFLEN);
	~WinSocks();

	//Non-blocking call to Listen for incoming messages - return message when received
	std::string Listen(void);

	//set data to be sent to client - not intended to send massive amounts of data, e.g. shouldn't use it to send the entire mesh data
	//SetSendData just loads values in dataParams std::vector whilst the SendDataParams actually sends data
	void SetSendData(const std::vector<std::string>& newdataParams);
	void SetSendData(std::vector<std::string>&& newdataParams);

	void SendDataParams(void);

	//if something went wrong with a received command then set the first data field as the '!' symbol - this will let the client know something is not right
	void SignalClientError(std::string message);
};

//-----------------------------------------------------------------------------------------------------------------------------------------------

//SERVER------------------------------------------------------------------------

inline WinSocks::WinSocks(int buffer_size) 
{
	recvbuflen = buffer_size;
	recvbuf = new char[recvbuflen];

	MakeListenSocket();
}

inline WinSocks::~WinSocks() {

	if (clientConnected) {

		// shutdown the connection since we're done
		int iResult = shutdown(ClientSocket, SD_SEND);
		if (iResult == SOCKET_ERROR) {

			//error
			closesocket(ClientSocket);
			WSACleanup();
		}

		// cleanup
		closesocket(ClientSocket);
	}

	if (listenSocketActive) {

		closesocket(ListenSocket);
	}

	WSACleanup();

	delete recvbuf;
}

inline void WinSocks::MakeListenSocket(void) {

	clientConnected = false;
	listenSocketActive = false;

	ListenSocket = INVALID_SOCKET;
	ClientSocket = INVALID_SOCKET;

	recvbuflen = DEFAULT_BUFLEN;

	// Initialize Winsock
	WSADATA wsaData;
	int iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
	if (iResult != 0) {

		//error
		return;
	}

	struct addrinfo *result = nullptr;
	struct addrinfo hints;

	ZeroMemory(&hints, sizeof(hints));
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_protocol = IPPROTO_TCP;
	hints.ai_flags = AI_PASSIVE;

	// Resolve the server address and port
	iResult = getaddrinfo(nullptr, DEFAULT_PORT, &hints, &result);
	if (iResult != 0) {

		//error
		WSACleanup();
		return;
	}

	// Create a SOCKET for connecting to server
	ListenSocket = socket(result->ai_family, result->ai_socktype, result->ai_protocol);
	if (ListenSocket == INVALID_SOCKET) {

		//error
		freeaddrinfo(result);
		WSACleanup();
		return;
	}

	// Setup the TCP listening socket (call ::bind, which is defined in Winsock2.h, not bind, which is defined in std)
	iResult = ::bind(ListenSocket, result->ai_addr, (int)result->ai_addrlen);
	if (iResult == SOCKET_ERROR) {

		//error
		freeaddrinfo(result);
		closesocket(ListenSocket);
		WSACleanup();
		return;
	}

	freeaddrinfo(result);

	iResult = listen(ListenSocket, SOMAXCONN);
	if (iResult == SOCKET_ERROR) {

		//error
		closesocket(ListenSocket);
		WSACleanup();
		return;
	}

	//Set ListenSocket as non-blocking (iMode = 1 for non-blocking, iMode = 0 for blocking)
	u_long iMode = 1;
	ioctlsocket(ListenSocket, FIONBIO, &iMode);

	listenSocketActive = true;
}

inline std::string WinSocks::Listen(void) {

	if (!listenSocketActive) { Sleep(RECVSLEEPMS); MakeListenSocket(); return ""; }

	// Accept a client socket if one not available already. Asynchronous since ListenSocket was made non-blocking.
	if (!clientConnected) {

		ClientSocket = accept(ListenSocket, nullptr, nullptr);

		if (ClientSocket != INVALID_SOCKET) {

			//Set ClientSocket as non-blocking (iMode = 1 for non-blocking, iMode = 0 for blocking)
			u_long iMode = 1;
			ioctlsocket(ClientSocket, FIONBIO, &iMode);

			//mark new client connection established
			clientConnected = true;
		}
		else {

			Sleep(RECVSLEEPMS);
			return "";
		}
	}

	//client is connected : listen for message
	int iResult = recv(ClientSocket, recvbuf, recvbuflen, 0);

	if (iResult == 0) {

		//Client connection was closed
		closesocket(ClientSocket);
		clientConnected = false;
		return "";
	}

	//WSAEWOULDBLOCK is issued for non-blocking sockets so ignore it, but respond to all other errors
	int nError = WSAGetLastError();
	if (nError && nError != WSAEWOULDBLOCK) {

		//error
		// Shutdown our socket
		shutdown(ClientSocket, SD_SEND);
		closesocket(ClientSocket);
		clientConnected = false;
		return "";
	}

	//all good, return the received message
	if (iResult > 0) {

		//received message, return it to be processed
		last_message = std::string(recvbuf).substr(0, iResult);

		return last_message;
	}

	//sleep before calling the non-blocking recv function (client socket is not blocking)
	Sleep(RECVSLEEPMS);

	return "";
}

inline void WinSocks::SetSendData(const std::vector<std::string>& newdataParams) {

	dataParams = newdataParams;
}

inline void WinSocks::SetSendData(std::vector<std::string>&& newdataParams) {

	dataParams = move(newdataParams);
}

inline void WinSocks::SendDataParams(void) {

	if (clientConnected) {

		//form message std::string : tab value tab value tab value etc. If no values then just a tab
		std::stringstream ss;

		if (!dataParams.size()) ss << "\t";

		for (int n = 0; n < dataParams.size(); n++) {

			ss << "\t" << dataParams[n];
		}

		//message client
		int iSendResult = send(ClientSocket, ss.str().c_str(), (int)ss.str().length(), 0);
		if (iSendResult == SOCKET_ERROR) {
			//error
			return;
		}
	}

	//done so reset std::vector to zero size
	dataParams.resize(0);
}

inline void WinSocks::SignalClientError(std::string message) {

	dataParams.resize(0);
	dataParams.push_back("!");
	dataParams.push_back(message);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//CLIENT

class WinSocksClient {

	WSADATA wsaData;
	SOCKET ConnectSocket;
	struct addrinfo *result, *ptr, hints;

	char *recvbuf;
	int recvbuflen;

public:

	bool connected;

public:

	WinSocksClient(std::string serverAddress = "", std::string port = DEFAULT_PORT, int buffer_size = DEFAULT_BUFLEN);
	~WinSocksClient();

	//Send message to server then wait for response
	std::string SendMessage_GetResponse(std::string message);

	// Send simple message to server
	void SendSimpleMessage(std::string message);
};

inline WinSocksClient::WinSocksClient(std::string serverAddress, std::string port, int buffer_size)
{
	//for local host use serverAddress = "localhost", or just leave it blank

	if (!serverAddress.length()) serverAddress = std::string("localhost");

	connected = false;

	ConnectSocket = INVALID_SOCKET;
	result = NULL;
	ptr = NULL;

	recvbuflen = buffer_size;
	recvbuf = new char[recvbuflen];

	// Initialize Winsock
	int iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
	if (iResult != 0) {
		
		return;
	}

	ZeroMemory(&hints, sizeof(hints));
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_protocol = IPPROTO_TCP;

	// Resolve the server address and port
	iResult = getaddrinfo(serverAddress.c_str(), port.c_str(), &hints, &result);
	if (iResult != 0) {
		
		WSACleanup();
		return;
	}

	// Attempt to connect to an address until one succeeds
	for (ptr = result; ptr != NULL; ptr = ptr->ai_next) {

		// Create a SOCKET for connecting to server
		ConnectSocket = socket(ptr->ai_family, ptr->ai_socktype, ptr->ai_protocol);
		if (ConnectSocket == INVALID_SOCKET) {
			
			WSACleanup();
			return;
		}

		// Connect to server.
		iResult = connect(ConnectSocket, ptr->ai_addr, (int)ptr->ai_addrlen);
		if (iResult == SOCKET_ERROR) {

			closesocket(ConnectSocket);
			ConnectSocket = INVALID_SOCKET;
			continue;
		}
		break;
	}

	freeaddrinfo(result);

	if (ConnectSocket == INVALID_SOCKET) {
		
		WSACleanup();
		return;
	}

	//now connected to server
	connected = true;
}

inline WinSocksClient::~WinSocksClient()
{
	// shutdown the connection since no more data will be sent
	int iResult = shutdown(ConnectSocket, SD_SEND);
	if (iResult == SOCKET_ERROR) {
		
		closesocket(ConnectSocket);
		WSACleanup();
	}

	//cleanup
	closesocket(ConnectSocket);
	WSACleanup();

	delete recvbuf;
}

inline std::string WinSocksClient::SendMessage_GetResponse(std::string message)
{
	if (!connected) return "";

	//tell server we're ready to receive message
	SendSimpleMessage(std::string(message));

	//now wait for message to be received
	int iResult = recv(ConnectSocket, recvbuf, recvbuflen, 0);
	if (iResult > 0) {

		std::string receivedString = std::string(recvbuf).substr(0, iResult);

		return receivedString;
	}
	else if (iResult == 0) {
		
		connected = false;
		return "";
	}
	else {

		connected = false;
		return "";
	}

	connected = false;
	return "";
}

inline void WinSocksClient::SendSimpleMessage(std::string message)
{
	// Send simple message

	int iResult = send(ConnectSocket, message.c_str(), message.length(), 0);
	if (iResult == SOCKET_ERROR) {

		closesocket(ConnectSocket);
		WSACleanup();
		connected = false;
		return;
	}
}