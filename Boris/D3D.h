//Defines low level Direct 3D methods

#pragma once

#include "CompileFlags.h"
#if OPERATING_SYSTEM == OS_WIN

#include <Wincodec.h>
#include <Windows.h>
#include <math.h>
#include <D3D11_1.h>
#include <d3dcompiler.h>
#include <directxcolors.h>
#include <DirectXMath.h>
#include <mutex>
#include <omp.h>
#include <mfapi.h>
#include <mfidl.h>
#include <Mfreadwrite.h>
#include <mferror.h>
#include <comdef.h>
#include <d2d1_1.h>
#include <dwrite_1.h>
#include <string>
#include <vector>

#include "BorisLib.h"

#include "ScreenGrab.h"

#pragma comment(lib, "mfreadwrite")
#pragma comment(lib, "mfplat")
#pragma comment(lib, "mfuuid")

#pragma comment (lib, "D3D11.lib")
#pragma comment (lib, "d2d1.lib")
#pragma comment (lib, "dwrite.lib")
#pragma comment (lib, "d3dcompiler.lib")

#define FARCLIPPINGPLANEDISTANCE 10000.0F

//pre-set video bit rates for encoding at the reference fps - scale the bit-rate used with fps to keep the same quality
#define VIDEO_BIT_RATE_VLOW	1000000
#define VIDEO_BIT_RATE_LOW	2500000
#define VIDEO_BIT_RATE_MEDIUM 4500000
#define VIDEO_BIT_RATE_HIGH 7000000
#define VIDEO_BIT_RATE_VHIGH 9500000
#define VIDEO_BIT_RATE_ULTRA 14000000
#define VIDEO_BIT_RATE_REFFPS	33

//Camera default settings
#define CAMERADISTANCE	300.0f
#define CAMERAFOV	30.0f

using namespace DirectX;


//enum defining cell display objects. These are custom objects so better not to define them here, but in the file implementing the derived client class (whose base class is D3D).
enum CDO_;

template<class Interface>
inline void SafeRelease(Interface **ppInterfaceToRelease)
{
	if (*ppInterfaceToRelease != nullptr) {

		(*ppInterfaceToRelease)->Release();
		(*ppInterfaceToRelease) = nullptr;
	}
}

struct SimpleVertex {
    
	XMFLOAT3 Pos;
    XMFLOAT4 Color;

	SimpleVertex(XMFLOAT3 _Pos, XMFLOAT4 _Color) {

		Pos = _Pos;
		Color = _Color;
	}

	SimpleVertex(void) {}
};

//Object constant buffer
struct ConstantBuffer {
	XMMATRIX mWorld;
	XMMATRIX mView;
	XMMATRIX mProjection;
	XMFLOAT4 vOutputColor;
};

////////////////////////////////////////////////////////////////////////////////////////////////

//Wrapper for all the stuff used to render an object. These are stored in ObjectBufferCollection for convenience.
//Make new ObjectBuffer object using ObjectBufferCollection::NewObjectBuffer method, built from vertices and indices, vertex and pixel shaders.
//Before drawing this object you must set context with ObjectBuffer::SetContext method - see D3D::DrawCBObjectBatch method. ObjectBufferCollection::SetContext should be used, as this will only set context for the required obbject if not already set
struct ObjectBuffer {

	ID3D11Buffer* pCB;
	
	//iCount has the number of indices to draw - use it with ID3D11DeviceContext::DrawIndexed method.
	UINT iCount;

	ID3D11VertexShader*     pVertexShader;
	ID3D11PixelShader*      pPixelShader;

	ID3D11InputLayout*		pVertexLayout;

	ID3D11Buffer*			pVertexBuffer;
	UINT stride, offset;

	ID3D11Buffer*			pIndexBuffer;

	ObjectBuffer(void) {

		iCount = 0;
		pCB = nullptr;

		pVertexShader = nullptr;
		pPixelShader = nullptr;

		pVertexLayout = nullptr;
		
		pVertexBuffer = nullptr;
		stride = 0; offset = 0;

		pIndexBuffer = nullptr;
	}

	~ObjectBuffer() { Release(); }

	void Release(void) {

		if(pCB) pCB->Release();
		pCB = nullptr;
		iCount = 0;

		SafeRelease(&pVertexShader);
		pVertexShader = nullptr;
		SafeRelease(&pPixelShader);
		pPixelShader = nullptr;

		SafeRelease(&pVertexLayout);
		pVertexLayout = nullptr;

		SafeRelease(&pVertexBuffer);
		pVertexBuffer = nullptr;
		stride = 0; offset = 0;

		SafeRelease(&pIndexBuffer);
		pIndexBuffer = nullptr;
	}

	void SetContext(ID3D11DeviceContext* pImmediateContext) {

		pImmediateContext->IASetInputLayout(pVertexLayout);

		pImmediateContext->IASetVertexBuffers(0, 1, &pVertexBuffer, &stride, &offset);

		pImmediateContext->IASetIndexBuffer(pIndexBuffer, DXGI_FORMAT_R16_UINT, 0);

		pImmediateContext->VSSetShader(pVertexShader, nullptr, 0);
		pImmediateContext->VSSetConstantBuffers(0, 1, &pCB);
		pImmediateContext->PSSetShader(pPixelShader, nullptr, 0);
		pImmediateContext->PSSetConstantBuffers(0, 1, &pCB);
	}
};

struct ObjectBufferCollection {

	std::vector<ObjectBuffer*> objBuf;

	CDO_ activeOB;

	ObjectBufferCollection(CDO_ activeOB) { this->activeOB = activeOB; }

	HRESULT NewObjectBuffer(ID3D11Device* pd3dDevice, ID3D11DeviceContext* pImmediateContext, SimpleVertex vertices[], UINT vSize, WORD indices[], UINT iSize, std::string VSFile, std::string PSFile);

	//Load shaders from a .fx file and compile it to a .cso file
	HRESULT CompileShaderFromFile(std::string fileName, LPCSTR szEntryPoint, LPCSTR szShaderModel, ID3DBlob** ppBlobOut);

	void SetContext(ID3D11DeviceContext* pImmediateContext, CDO_ objColSelector) {

		if(activeOB != objColSelector) {

			activeOB = objColSelector;

			objBuf[objColSelector]->SetContext(pImmediateContext);
		}
	}

	void Release(void) {

		for(int n = 0; n < objBuf.size(); n++) {

			objBuf[n]->Release();
		}
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////

//used to store transformations for a constant buffer object - rotation, scaling, translation and color. Make a vector with this to draw a batch of identical objects.
struct CBObjectTransform {

	XMMATRIX Rotat;
	XMMATRIX Scale;
	XMMATRIX Trans;
	XMFLOAT4 Color;

	//maximum set transparency (alpha value)
	double maxTransparency = 1.0;

	//the "value" of the quantity for which this transform was set. e.g. if a scalar physical quantity then this is its value, if a vector then this is its magnitude.
	double value;

	//the position of object (in world view, logical units) on which this Transform is used : this will be the translation value used to set Trans
	DBL3 position;

	//set this to true to skip rendering (e.g. empty or obscured cell)
	bool skip_render = false;

	CBObjectTransform(void) 
	{
		Rotat = XMMatrixIdentity();
		Scale = XMMatrixIdentity() * XMMatrixScaling(0, 0, 0);
		Trans = XMMatrixIdentity();
		Color = XMFLOAT4(0, 0, 0, 1);
		value = 0.0;
		position = DBL3();
	}

	CBObjectTransform(XMMATRIX _Rotat, XMMATRIX _Scale, XMMATRIX _Trans, XMFLOAT4 _Color, double _maxTransparency, DBL3 position_, double _value) 
	{
		Rotat = _Rotat;
		Scale = _Scale;
		Trans = _Trans;
		Color = _Color;
		maxTransparency = _maxTransparency;
		value = _value;
		position = position_;
	}

	CBObjectTransform(const CBObjectTransform& copyThis) { *this = copyThis; }
	CBObjectTransform& operator=(const CBObjectTransform& copyThis)
	{
		Rotat = copyThis.Rotat;
		Scale = copyThis.Scale;
		Trans = copyThis.Trans;
		Color = copyThis.Color;
		maxTransparency = copyThis.maxTransparency;

		value = copyThis.value;
		position = copyThis.position;

		return *this;
	}

	inline XMMATRIX Transform(void) 
	{	 
		return Rotat * Scale * Trans;
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////

class D3D {

private:

	//world, viewspace and projection space matrices
	XMMATRIX g_World, g_View, g_Projection;

	D3D11_VIEWPORT viewPort;

	//camera settings. Eye: position of camera in world space, At: focal point of camera, Up: the up direction of camera view
	XMVECTOR Eye, At, Up;

	//pixels per unit height of 3D viewport, i.e. the screen height in pixels divided by the camera viewport height
	float pixperUnit;

	//Initialized : d3d can be used. d2dInitialized : d2d can also be used.
	bool Initialized, d2dInitialized;

	D3D_DRIVER_TYPE         driverType;
	D3D_FEATURE_LEVEL       featureLevel;

	//dxgi factory used to make things, among them the swap chain
	IDXGIFactory1*			pdxgiFactory;

	//swap chain : used to present back buffer to front buffer (screen)
	IDXGISwapChain*         pSwapChain;
	//the back buffer 
	ID3D11RenderTargetView* pRenderTargetView;

	//depth buffer : enabled to allow depth rendering (i.e. farther away objects are always obscured by closer objects, irrespective of order of drawing)
	ID3D11DepthStencilView* pDepthStencilView;

protected:

	//camera position, distance from origin and fov in degrees
	float camX, camY, camZ, camDistance, fovDeg;

	//camera up direction in Cartesian coordinates
	float camUx, camUy, camUz;

	//shift g_View by this (pixel values) - translates Eye, At, Up. 0 values result in view being set to screen center
	float view_shiftX, view_shiftY, view_shiftZ;

	//dpi settings of screen
	float dpiX, dpiY;

	//d3d device: used to make things
	ID3D11Device*           pd3dDevice;

	//device context : used to render things (to RenderTargetView, a.k.a. back buffer) before being presented
	ID3D11DeviceContext*    pImmediateContext;

	//d2d factories and render target
	ID2D1RenderTarget*		pD2DRT;
	ID2D1Factory*			pD2DFactory;
	IDWriteFactory1*		pWriteFactory;

public:  //Public data

	//main window width and height
	UINT wndWidth, wndHeight;

	//Use a std::mutex to access drawing methods which could be accessed by multiple threads. 
	//Only use the std::mutex at the first point of entry : define a Refresh method which handles all calls to updating the screen (parts or all of it). The std::mutex is used there.
	//Also use the std::mutex at the start of methods which call for drawing to be done to a file, as these will also use the back buffer.
	//NOTE : this std::mutex may not be needed, if instead graphics routines are managed by a top level display object (calls eventually get here), in which case it should have its own std::mutex.
	std::mutex GraphicsMutex;

private: //Private methods

	//Clean up all allocated resources.
	void CleanupDevice(void);

	//Initialize d3d and d2d interop respectively
	HRESULT InitDevice(HWND hWnd);
	HRESULT SetD2DInterop(void);

	//calculate conversion constant between screen pixels and screen-projected 3D view logical units
	void CalculatePixPerUnit(void);

protected: //Protected methods

	//begin and end drawing using D2D render target. All D2D drawing MUST be enclosed between this pair.
	//The procedure used here is to have a Refresh method which starts and ends with the BeginD3DDraw() ... EndD3DDraw() pair. Any D2D methods used should individually have a BeginD2DDraw() ... EndD2DDraw() pair.
	inline void BeginD2DDraw(void) {pD2DRT->BeginDraw();}
	inline void EndD2DDraw(void) {pD2DRT->EndDraw();}

	//get std::string from HRESULT
	std::string GetErrorMessage(HRESULT hr);

	//Draw a batch of identical objects (objBuf with given pixel and vertex shaders) as indicated by the cbOTBatch vector - each object can be individually rotated, scaled, translated and colored
	void DrawCBObjectBatch(std::vector<CBObjectTransform> &cbOTBatch, ObjectBufferCollection &objCol, CDO_ objColSelector);
	//Draw a single constant buffer object
	void DrawCBObject(XMMATRIX Rotation, XMMATRIX Scaling, XMMATRIX Translation, XMFLOAT4 Color, ObjectBufferCollection &objCol, CDO_ objColSelector);

	//Draw a single object using linelist primitive topology
	void DrawFrameObject(XMMATRIX Rotation, XMMATRIX Scaling, XMMATRIX Translation, XMFLOAT4 Color, ObjectBufferCollection &objCol, CDO_ objColSelector);

	//save currently displaying image in capture_rect to png file
	bool SaveScreenToPNG(std::string fileName, D2D1_RECT_F capture_rect);

	//Load bitmap from file to given size
	HRESULT LoadScaledBitmap(std::string fileName, IWICFormatConverter **ppWICFormatConverter, D2D1_SIZE_U size);
	//Get dimensions in pixels of image in the given image file
	HRESULT GetImageDimensions(std::string fileName, UINT *pbmpWidth, UINT *pbmpHeight);

	//Video encoding methods : Initialize a sink writer and write frames to it. These methods must be directed by another method which starts the process, writes frames and ends the process.
	//Initialize a sink writer and also get a stream index for it. The output video will be in the given filename, with the given fps and frames will have the given dimensions.
	//quality determines the bit-rate used
	HRESULT InitializeSinkWriter(std::string fileName, IMFSinkWriter **ppWriter, DWORD *pStreamIndex, UINT32 VIDEO_WIDTH, UINT32 VIDEO_HEIGHT, UINT32 VIDEO_FPS, int quality);
	//Write a frame to the initialized sink writer (also need the stream index for it) by loading it from the given image file. The frame is loaded to the given dimensions (which must match the initialized sink writer).
	//The fps must also match that of the initialized sink writer. The frame will be placed at the time position indicated by timeStamp: each frame lasts for VIDEO_FRAME_DURATION = 10 * 1000 * 1000 / VIDEO_FPS ticks, so increment timeStamp by this count for every new frame.
	HRESULT WriteFrameFromFile(std::string fileName, IMFSinkWriter *pWriter, DWORD streamIndex, const LONGLONG& timeStamp, UINT32 VIDEO_WIDTH, UINT32 VIDEO_HEIGHT, UINT32 VIDEO_FPS);

public:  //Public methods

	D3D(HWND hWnd);
	virtual ~D3D();

	//Any rendering should be done between a pair of BeginD3DDraw() ... EndD3DDraw(), including drawing using the D2D render target (pD2DRT), but the latter must also use the D2D pairing methods.
	//clears screen and other buffers as necessary - we just clear the depth buffer here.
	void BeginD3DDraw(void) { pImmediateContext->ClearDepthStencilView(pDepthStencilView, D3D11_CLEAR_DEPTH, 1.0f, 0); }
	void ClearScreen(void) { pImmediateContext->ClearRenderTargetView(pRenderTargetView, Colors::White); }
	//presents back buffer to front buffer (screen)
	void EndD3DDraw(void) { pSwapChain->Present( 0, 0 ); }

	//call ReInitialize whenever window dimensions have changed. It will automatically clean and remake all resources held in this class.
	void ReInitialize(HWND hWnd);

	/////////////////////Methods used to control camera/////////////////////////
	
	void SetCameraPosition(float camX, float camY, float camZ);
	void SetCameraUp(float camUx, float camUy, float camUz);
	void SetCameraFOV(float fovDeg);

	//Combination of SetCameraPosition, SetCameraUp and  SetCameraFOV
	void SetCamera(float camX, float camY, float camZ, float camUx, float camUy, float camUz, float fovDeg);
	
	//set to new camera distance
	void SetCameraDistance(float camDistance);
	
	//Set the shift from the screen centre of the on-screen origin of 3D system by given pixel distances
	void Shift3DOriginPixelPosition(float dX, float dY, float dZ = 0);
	
	//Set the on-screen origin of 3D system at given screen pixel position
	void Set3DOriginPixelPosition(float X, float Y, float Z = 0);

	//set camera to default settings
	void ResetCamera(void);

	//Make a video file from an input image sequence at the given fps. The input images are scaled in dimensions by the frameScaler value (best quality frameScaler = 1, worst quality capped at frameScaler = 0.25)
	bool MakeVideoFromFileSequence(std::string outputFile, std::vector<std::string> &inputFiles, UINT32 VIDEO_FPS, float frameScaler, int quality);

	////////////////Getters

	float GetPixPerUnit(void) { return pixperUnit; }

	FLT3 GetCameraPosition(void) { return FLT3(camX, camY, camZ); }
	FLT3 GetCameraUp(void) { return FLT3(camUx, camUy, camUz); }
	float GetCameraFOV(void) { return fovDeg; }
	float GetCameraDistance(void) { return camDistance; }

	//from mouse coordinates on screen return corresponding world position on far clipping plane (in logical units)
	DBL3 Pick_FarPlane_Point(INT2 mouse);

	////////////////Work with image files

	//from image file (.png) extract a raw bitmap (BYTE array with 4 BYTE-sized entries as B-G-R-A for each pixel) to specified pixel size (so image file is rescaled to specified size)
	bool GetBitmapFromImage(std::string bitmapFile, std::vector<unsigned char>& bitmap, INT2 pixels_size);

};

#endif