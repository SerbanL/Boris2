#include "stdafx.h"
#include "WinSocks.h"

//SERVER------------------------------------------------------------------------

WinSocks::WinSocks(void) {

	MakeListenSocket();
}

WinSocks::~WinSocks() {

	if(clientConnected) {

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

	if(listenSocketActive) {

		closesocket(ListenSocket);
	}

	WSACleanup();
}

void WinSocks::MakeListenSocket(void) {

	clientConnected = false;
	listenSocketActive = false;

	ListenSocket = INVALID_SOCKET;
    ClientSocket = INVALID_SOCKET;

    recvbuflen = DEFAULT_BUFLEN;

	// Initialize Winsock
	WSADATA wsaData;
    int iResult = WSAStartup(MAKEWORD(2,2), &wsaData);
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
    if ( iResult != 0 ) {
        
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

    // Setup the TCP listening socket
    iResult = bind(ListenSocket, result->ai_addr, (int)result->ai_addrlen);
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
	ioctlsocket(ListenSocket,FIONBIO,&iMode);

	listenSocketActive = true;
}

string WinSocks::Listen(void) {

	if(!listenSocketActive) { Sleep(RECVSLEEPMS); MakeListenSocket(); return ""; }

	// Accept a client socket if one not available already. Asynchronous since ListenSocket was made non-blocking.
	if(!clientConnected) {

		ClientSocket = accept(ListenSocket, nullptr, nullptr);

		if(ClientSocket != INVALID_SOCKET) {

			//Set ClientSocket as non-blocking (iMode = 1 for non-blocking, iMode = 0 for blocking)
			u_long iMode = 1;
			ioctlsocket(ClientSocket,FIONBIO,&iMode);

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

	if(iResult == 0) {

		//Client connection was closed
		closesocket(ClientSocket);
		clientConnected = false;
        return "";
	}

	//WSAEWOULDBLOCK is issued for non-blocking sockets so ignore it, but respond to all other errors
	int nError = WSAGetLastError();
	if(nError && nError != WSAEWOULDBLOCK) {

		//error
		// Shutdown our socket
		shutdown(ClientSocket,SD_SEND);
		closesocket(ClientSocket);
		clientConnected = false;
        return "";
	}

	//all good, return the received message
	if (iResult > 0) {

		//received message, return it to be processed
		string receivedString = string(recvbuf).substr(0, iResult);

		return receivedString;
    }

	//sleep before calling the non-blocking recv function (client socket is not blocking)
	Sleep(RECVSLEEPMS);

	return "";
}

void WinSocks::SetSendData(const vector<string>& newdataParams) {

	dataParams = newdataParams;
}

void WinSocks::SetSendData(vector<string>&& newdataParams) {
	
	dataParams = move( newdataParams );
}

void WinSocks::SendDataParams(void) {
	
	if(clientConnected) {
		
		//form message string : tab value tab value tab value etc. If no values then just a tab
		stringstream ss;

		if(!dataParams.size()) ss << "\t";
		
		for(int n = 0; n < dataParams.size(); n++) {

			ss << "\t" << dataParams[n];
		}

		//message client
		int iSendResult = send( ClientSocket, ss.str().c_str(), (int)ss.str().length(), 0 );
		if(iSendResult == SOCKET_ERROR) {
			//error
			return;
		}
	}
	
	//done so reset vector to zero size
	dataParams.resize(0);
}

void WinSocks::SignalClientError(string message) {

	dataParams.resize(0);
	dataParams.push_back("!");
	dataParams.push_back(message);
}