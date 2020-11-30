#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

//default text color options
//See: https://msdn.microsoft.com/en-us/library/windows/desktop/dd370907(v=vs.85).aspx

#define MESSAGECOLOR D2D1::ColorF(D2D1::ColorF::Green)  //rgba = (0, 0.5, 0, 1)
#define ERRORCOLOR D2D1::ColorF(D2D1::ColorF::Red)		//rgba = (1, 0, 0, 1)
#define WARNINGCOLOR D2D1::ColorF(D2D1::ColorF::Yellow)	//rgba = (1, 1, 0, 1)
#define LISTINGCOLOR D2D1::ColorF(D2D1::ColorF::White)  //rgba = (1, 1, 1, 1)
#define USERCOLOR D2D1::ColorF(D2D1::ColorF::Yellow)	//rgba = (1, 1, 0, 1)
#define BGRNDTEXTCOLOR D2D1::ColorF(0,0)				//rgba = (0, 0, 0, 0)
#define PROMPTCOLOR D2D1::ColorF(13882323, 0.7)			//Light Gray: rgba = (211/255, 211/255, 211/255, 1). Thus single rgba value as UINT32 : 13882323 = 256^2 * 211 + 256 * 211 + 211
#define POPUPCOLOR D2D1::ColorF(11119017, 0.7)			//Gray: rgba = (255/255, 255/255, 255/255, 1). Thus single rgba value as UINT32 : 13882323 = 256^2 * 255 + 256 * 255 + 255
#define HOVERINFOCOLOR D2D1::ColorF(D2D1::ColorF::Blue)
#define ONCOLOR D2D1::ColorF(D2D1::ColorF::Green)
#define ALTONCOLOR D2D1::ColorF(D2D1::ColorF::Orange)
#define OFFCOLOR D2D1::ColorF(D2D1::ColorF::Red)
#define UNAVAILABLECOLOR D2D1::ColorF(D2D1::ColorF::Gray)
#define UNAVAILABLECOLOR2 D2D1::ColorF(D2D1::ColorF::DarkGray)
#define RESIZEFRAMECOLOR D2D1::ColorF(D2D1::ColorF::Orange)
#define FRAMECOLOR D2D1::ColorF(D2D1::ColorF::Yellow)

#else

//In text console mode these are meaningless but need to keep rest of the code happy (and clean)

#define MESSAGECOLOR	0
#define ERRORCOLOR	0
#define LISTINGCOLOR	0
#define USERCOLOR	0
#define BGRNDTEXTCOLOR	0
#define PROMPTCOLOR	0
#define POPUPCOLOR	0
#define HOVERINFOCOLOR	0
#define ONCOLOR	0
#define OFFCOLOR	0
#define UNAVAILABLECOLOR	0
#define UNAVAILABLECOLOR2	0
#define RESIZEFRAMECOLOR	0

#endif
