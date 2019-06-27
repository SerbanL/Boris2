#include "stdafx.h"
#include "BorisWinapi.h"

#include "CompileFlags.h"
#if GRAPHICS == 1 && OPERATING_SYSTEM == OS_WIN

#define MAX_LOADSTRING 2000

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Program master object - WndProc dispatches messages to this.
Simulation *pSim = nullptr;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Global Variables:
HINSTANCE hInst;								// current instance
TCHAR szTitle[MAX_LOADSTRING];	 				// The title bar text
TCHAR szWindowClass[MAX_LOADSTRING];			// the main window class name

// Forward declarations of functions included in this code module:
ATOM				MyRegisterClass(HINSTANCE hInstance);
BOOL				InitInstance(HINSTANCE, int);
LRESULT CALLBACK	WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK	About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY _tWinMain(_In_ HINSTANCE hInstance,
					   _In_opt_ HINSTANCE hPrevInstance,
					   _In_ LPTSTR    lpCmdLine,
					   _In_ int       nCmdShow)
{
	UNREFERENCED_PARAMETER(hPrevInstance);
	UNREFERENCED_PARAMETER(lpCmdLine);

 	// TODO: Place code here.
	MSG msg;

	// Initialize global strings
	LoadString(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
	LoadString(hInstance, IDC_BORISWINAPI, szWindowClass, MAX_LOADSTRING);
	
	//overwrite the title bar

	//Bizarre bug here : if you convert to (float) instead of (double) the program crashes!!!
	//If however you declare a stringstream variable before (just declare, you don't have to use it!) then program doesn't crash!
	//Gets better : it used to work with (float) but for an unknown reason it suddenly started crashing and I don't know what's changed.
	//Without a doubt the strangest bug I have ever seen - but as usual it will make sense when investigated properly, but I suspect it won't be easy. TO DO if you have time.
	string sTitle = string("Boris v") + ToString((double)Program_Version / 100);
	copy(sTitle.begin(), sTitle.end(), szTitle);
	szTitle[sTitle.size()] = 0;

	MyRegisterClass(hInstance);

	// Perform application initialization: (use SW_MAXIMIZE to have the base window maximized)
	if (!InitInstance (hInstance, nCmdShow))
	{
		return FALSE;
	}

	HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_BORISWINAPI));
	
	//Windows message loop
	while(GetMessage(&msg, nullptr, 0, 0)) {

		if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg)) {

			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}
	
	return (int) msg.wParam;
}

//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
	WNDCLASSEX wcex;

	wcex.cbSize = sizeof(WNDCLASSEX);

	wcex.style			= CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc	= WndProc;
	wcex.cbClsExtra		= 0;
	wcex.cbWndExtra		= 0;
	wcex.hInstance		= hInstance;
	wcex.hIcon			= LoadIcon(hInstance, MAKEINTRESOURCE(IDI_BORISWINAPI));
	wcex.hCursor		= LoadCursor(nullptr, IDC_ARROW);
	wcex.hbrBackground	= (HBRUSH)(COLOR_WINDOW+1);
	//take out the following so we don't have a menu (default menu includes File and About options).
	//wcex.lpszMenuName	= MAKEINTRESOURCE(IDC_BORISWINAPI);
	wcex.lpszClassName	= szWindowClass;
	wcex.hIconSm		= LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

	return RegisterClassEx(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   HWND hWnd;

   hInst = hInstance; // Store instance handle in our global variable

   hWnd = CreateWindow(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW | WS_EX_ACCEPTFILES,
	   CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, SW_SHOWMAXIMIZED);
   UpdateWindow(hWnd);
   
   //Stop the "Program has stopped working" pop-up error box after closing program. I really cannot figure out why it appears sometimes!!!
   //!!! NOTE !!! If this is not enabled and the program closes with the "Program has stopped working" error, then the exe file is put on some windows watch list
   //(some kind of hook inserted in the running code to check function calls?), which will slow down execution dramatically. 
   //It seems this is done by the DiagTrack service (Diagnostics Tracking Service) - stop DPS service (Diagnostic Policy Service) then DiagTrack.
   SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);

   //Set normal arrow cursor - we are answering WM_SETCURSOR message and doing nothing there so need to set cursor separately
   SetCursor(LoadCursor(nullptr, IDC_ARROW));
   
   //Set to accept dropped files - file drop handled by WM_DROPFILES
   DragAcceptFiles(hWnd, true);

   //Allow WM_DROPFILES through UIPI filter for Windows 10 : if running in Administrator mode this message will be filtered since drag and drop from an unelevated source to an elevated recipient is blocked by default
   //This is what you have to do to force Windows 10 to allow drag and drop (marvellous isn't it!):
   ChangeWindowMessageFilterEx(hWnd, WM_DROPFILES, MSGFLT_ALLOW, nullptr);
   ChangeWindowMessageFilterEx(hWnd, WM_COPYDATA, MSGFLT_ALLOW, nullptr);
   unsigned int WM_COPYGLOBALDATA = 0x0049;
   ChangeWindowMessageFilterEx(hWnd, WM_COPYGLOBALDATA, MSGFLT_ALLOW, nullptr);

   //Instantiate Simulation object and start main simulation loop thread (non-blocking)
   pSim = new Simulation(hWnd, Program_Version);
   
   return TRUE;
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	PAINTSTRUCT ps;
	HDC hdc = 0;
	
	switch (message) {
	
	case WM_PAINT:
		hdc = BeginPaint(hWnd, &ps);
		if(pSim) pSim->RefreshScreen();
		EndPaint(hWnd, &ps);
		break;

	case WM_CLOSE:
		{
			//clean up and exit
			if(pSim) delete pSim;
			pSim = nullptr;

			//Destroy window : WM_DESTROY will be issued
			DestroyWindow(hWnd);
		}
		break;

	case WM_DESTROY:
		//Quit program : WM_QUIT will be issued
		PostQuitMessage(0);
		break;
		
	//keyboard key down - treat some special character codes here
	case WM_KEYDOWN: 
		{
			
			//intercept CTRL+ something combinations - NOTE (!!!), if that something is an ascii character, it will also raise the WM_CHAR case, so you must check there that CTRL is not on
			if(GetKeyState(VK_CONTROL) & 0x8000) {
				
				switch(wParam) {

					case 0x56: //Ctrl+V
						{
							if(pSim) pSim->NewMessage(AC_CTRL_V, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), GetClipboardText());
						}
						break;

					case 0x43: //Ctrl+C
						{
							if(pSim) pSim->NewMessage(AC_CTRL_C, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						}
						break;

					default:
						break;
				}
			}
			else {
				
				switch (wParam) { 
			
				case VK_LEFT:   // LEFT ARROW
					if(pSim) {

						if(!(GetKeyState(VK_SHIFT) & 0x8000)) pSim->NewMessage(AC_KBDLEFT, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						else pSim->NewMessage(AC_KBDSHIFTLEFT, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					}
					break;
			
				case VK_RIGHT:  // RIGHT ARROW
					if(pSim) {

						if(!(GetKeyState(VK_SHIFT) & 0x8000)) pSim->NewMessage(AC_KBDRIGHT, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						else pSim->NewMessage(AC_KBDSHIFTRIGHT, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					}
					break; 

				case VK_UP:     // UP ARROW
					if(pSim) pSim->NewMessage(AC_KBDUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;
			
				case VK_DOWN:   // DOWN ARROW 
					if(pSim) pSim->NewMessage(AC_KBDDN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;
			
				case VK_HOME:   // HOME
					if(pSim) {

						if(!(GetKeyState(VK_SHIFT) & 0x8000)) pSim->NewMessage(AC_KBDHOME, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						else pSim->NewMessage(AC_KBDSHIFTHOME, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					}
					break;
			
				case VK_END:    // END
					if(pSim) {

						if(!(GetKeyState(VK_SHIFT) & 0x8000)) pSim->NewMessage(AC_KBDEND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						else pSim->NewMessage(AC_KBDSHIFTEND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					}
					break;
			
				case VK_INSERT: // INS
					if(pSim) pSim->NewMessage(AC_KBDINS, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;
			
				case VK_DELETE: // DEL
					if(pSim) pSim->NewMessage(AC_KBDDEL, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_PRIOR:	//PAGE UP
					if(pSim) pSim->NewMessage(AC_KBDPGUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_NEXT:	//PAGE DOWN
					if(pSim) pSim->NewMessage(AC_KBDPGDN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_ESCAPE:	//ESC
					if(pSim) pSim->NewMessage(AC_KBDESC, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_BACK:	//BACKSPACE
					if(pSim) pSim->NewMessage(AC_KBDBACK, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_RETURN: //ENTER KEY
					if(pSim) pSim->NewMessage(AC_KBDENTER, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				default:
					break;
				}
			}
		}
		break;
			
	//treat ascii character codes
	case WM_CHAR:
		{

			//if an ascii character key is pressed, first make sure we don't have CTRL key on as CTRL+something combinations are handled elsewhere. If you don't check here, the something key will also be handled here - not good!
			if(!(GetKeyState(VK_CONTROL) & 0x8000)) {

				switch (wParam) {
			
				case 0x08: // Process a backspace. 
					break;
			
				case 0x0A: // Process a linefeed. 
					break;
			
				case 0x1B: // Process an escape. 
					break; 
			
				case 0x09: // Process a tab. 
					break; 
			
				case 0x0D: // Process a carriage return.
					break;
			
				default: // Process displayable characters. 			
					if(pSim) pSim->NewMessage(AC_KEYBOARD, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString((char)wParam));
					break; 
				}
			}
		}
		break;
			
	case WM_LBUTTONDBLCLK: //mouse left button double click
		//NOTE: Need CS_DBLCLKS style set when creating Instance, for it to issue WM_LBUTTONDBLCLK.
		if(pSim) pSim->NewMessage(AC_DOUBLECLICK, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_LBUTTONDOWN: //left mouse button down
		if (pSim) {
			
			if ((GetKeyState(VK_SHIFT) & 0x8000) == 0x8000) pSim->NewMessage(AC_SHIFT_MOUSELEFTDOWN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
			else pSim->NewMessage(AC_MOUSELEFTDOWN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		}
		break;

	case WM_LBUTTONUP:  //left mouse button up
		if(pSim) pSim->NewMessage(AC_MOUSELEFTUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_MBUTTONDOWN:	//middle mouse button down
		if(pSim) pSim->NewMessage(AC_MOUSEMIDDLEDOWN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_MBUTTONUP:	//middle mouse button up
		if(pSim) pSim->NewMessage(AC_MOUSEMIDDLEUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;
		
	case WM_RBUTTONDOWN: //right mouse button down
		if(pSim) pSim->NewMessage(AC_MOUSERIGHTDOWN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_RBUTTONUP: //right mouse button up
		if(pSim) pSim->NewMessage(AC_MOUSERIGHTUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_MOUSEMOVE:
		if(pSim) pSim->NewMessage(AC_MOUSEMOVE, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_MOUSEWHEEL: //mouse wheel used
		if(pSim) pSim->NewMessage(AC_MOUSEWHEEL, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(GET_WHEEL_DELTA_WPARAM(wParam)/120));
		break;
		
	case WM_SETCURSOR:
		//answering this message to stop cursor from reverting to the default IDC_ARROW setting
		break;
		
	case WM_DROPFILES:
		{
			TCHAR lpszFile[MAX_LOADSTRING] = { 0 };
			UINT uFile = 0;
			HDROP hDrop = (HDROP)wParam;

			uFile = DragQueryFile(hDrop, 0xFFFFFFFF, NULL, NULL);
			if (uFile != 1) {
			
				//need only one dropped file
				DragFinish(hDrop);
				break;
			}

			lpszFile[0] = '\0';
		
			if (DragQueryFile(hDrop, 0, lpszFile, MAX_LOADSTRING)) {
			
				if (pSim) {

					POINT dropPoint;
					DragQueryPoint(hDrop, &dropPoint);

					string fileName = WideStringtoString(lpszFile);

					pSim->NewMessage(AC_DROPFILES, INT2((int)dropPoint.x, (int)dropPoint.y), fileName);
				}
			}

			DragFinish(hDrop);
		}
		break;

	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
		break;
	}
	
	return 0;
}

#endif