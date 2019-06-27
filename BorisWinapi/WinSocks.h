#pragma once

//#include "BorisLib.h" -> cannot include this as it includes <functional>. There's a function std::bind in functional, and another function bind (used to setup a TCP listening socket) defined in winsock2.h. You could remove using namespace std here, but BorisLib will need to be modified also.

#include <windows.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>

#pragma comment (lib, "Ws2_32.lib")

#define DEFAULT_BUFLEN 512
#define DEFAULT_PORT "1542"

#define RECVSLEEPMS		100				//Used for non-blocking calls to recv : after each call sleep this number of ms to avoid hogging too much CPU time 
										//Calling recv too often, or using blocking sockets i.e. blocking recv, will put a slight load on one of the CPU cores which can really degrade parallelized computations, irrespective of which thread what

using namespace std;

//see https://msdn.microsoft.com/en-us/library/windows/desktop/ms737889(v=vs.85).aspx

//SERVER

class WinSocks {

    SOCKET ListenSocket;
    SOCKET ClientSocket;

	//receive buffer and length
    char recvbuf[DEFAULT_BUFLEN];
    int recvbuflen;

	//flags
	bool clientConnected;
	bool listenSocketActive;

	//commands processed by Simulation may return parameters. These are placed in this vector to be sent to Client in winsocks thread
	vector<string> dataParams;

private:

	//Make non-blocking listen socket on DEFAULT_PORT to accept incoming connections - done in the constructor
	void MakeListenSocket(void);

public:

	WinSocks(void);
	~WinSocks();

	//Non-blocking call to Listen for incoming messages - return message when received
	string Listen(void);

	//set data to be sent to client - not intended to send massive amounts of data, e.g. shouldn't use it to send the entire mesh data
	//SetSendData just loads values in dataParams vector whilst the SendDataParams actually sends data
	void SetSendData(const vector<string>& newdataParams);
	void SetSendData(vector<string>&& newdataParams);

	void SendDataParams(void);

	//if something went wrong with a received command then set the first data field as the '!' symbol - this will let the client know something is not right
	void SignalClientError(string message);
};