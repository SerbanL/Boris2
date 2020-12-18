#include "stdafx.h"
#include "D3D.h"
#if OPERATING_SYSTEM == OS_WIN

D3D::D3D(HWND hWnd) {

	//initializations
	driverType = D3D_DRIVER_TYPE_NULL;
	featureLevel = D3D_FEATURE_LEVEL_11_0;

	pd3dDevice = nullptr;
	pdxgiFactory = nullptr;
	pImmediateContext = nullptr;
	
	pSwapChain = nullptr;
	pRenderTargetView = nullptr;

	pDepthStencilView = nullptr;

	pD2DRT = nullptr;
	pD2DFactory = nullptr;
	pWriteFactory = nullptr;
	
	ReInitialize(hWnd);

	// Initialize the world matrix
	g_World = XMMatrixIdentity();

	// Initialize the view and projection matrices to default settings
	ResetCamera();
}

D3D::~D3D() {

	CleanupDevice();
}

void D3D::ReInitialize(HWND hWnd) {

	CleanupDevice();

	driverType = D3D_DRIVER_TYPE_NULL;
	featureLevel = D3D_FEATURE_LEVEL_11_0;

	if(SUCCEEDED(InitDevice(hWnd))) Initialized = true;
	if(SUCCEEDED(SetD2DInterop())) d2dInitialized = true;

	CalculatePixPerUnit();

	pD2DFactory->GetDesktopDpi(&dpiX, &dpiY);

    // Initialize the projection matrix
	g_Projection = XMMatrixPerspectiveFovRH( XMConvertToRadians(fovDeg), (FLOAT)wndWidth / wndHeight, 0.01f, FARCLIPPINGPLANEDISTANCE );
}

void D3D::CalculatePixPerUnit(void) {
	
	//pixels per unit height of 3D viewport, i.e. the screen height in pixels divided by the camera viewport height in logical units
	//camera viewport height = 2 * distance * tan (fov/2) : imagine the camera viewport as an isosceles triangle with fov being its tip angle.
	//use this to convert between camera logical units and number of pixels

	pixperUnit = wndHeight / (2*camDistance*tan(XMConvertToRadians(fovDeg/2)));
}

void D3D::SetCameraPosition(float camX, float camY, float camZ) {
	
	//Camera position
	this->camX = camX; this->camY = camY, this->camZ = camZ;
	camDistance = GetMagnitude(camX, camY, camZ);
	Eye = XMVectorSet(camX, camY, camZ, 0);
	g_View = XMMatrixLookAtRH(Eye, At, Up);

	CalculatePixPerUnit();
}

void D3D::SetCameraDistance(float camDistance) {

	if(IsZoN(camDistance) || IsGE(camDistance, FARCLIPPINGPLANEDISTANCE)) return;

	//Normalised camera position
	float camX_norm = camX/this->camDistance;
	float camY_norm = camY/this->camDistance;
	float camZ_norm = camZ/this->camDistance;
	
	//Set new camera position
	SetCameraPosition(camX_norm*camDistance, camY_norm*camDistance, camZ_norm*camDistance);
}

void D3D::SetCameraUp(float camUx, float camUy, float camUz) {

	//Camera up direction
	this->camUx = camUx; this->camUy = camUy; this->camUz = camUz;
	Up = XMVectorSet(camUx, camUy, camUz, 0);
	g_View = XMMatrixLookAtRH(Eye, At, Up);
}

void D3D::SetCameraFOV(float fovDeg) {
	
	this->fovDeg = fovDeg;
	g_Projection = XMMatrixPerspectiveFovRH( XMConvertToRadians(fovDeg), (FLOAT)wndWidth / wndHeight, 0.01f, FARCLIPPINGPLANEDISTANCE );

	CalculatePixPerUnit();
}

//Combination of SetCameraPosition, SetCameraUp and  SetCameraFOV
void D3D::SetCamera(float camX, float camY, float camZ, float camUx, float camUy, float camUz, float fovDeg) {

	//Camera position
	this->camX = camX; this->camY = camY, this->camZ = camZ;
	camDistance = GetMagnitude(camX, camY, camZ);
	Eye = XMVectorSet(camX, camY, camZ, 0);

	//Camera up direction
	this->camUx = camUx; this->camUy = camUy; this->camUz = camUz;
	Up = XMVectorSet(camUx, camUy, camUz, 0);
	g_View = XMMatrixLookAtRH(Eye, At, Up);

	//Camera FOV
	this->fovDeg = fovDeg;
	g_Projection = XMMatrixPerspectiveFovRH( XMConvertToRadians(fovDeg), (FLOAT)wndWidth / wndHeight, 0.01f, FARCLIPPINGPLANEDISTANCE );

	CalculatePixPerUnit();
}

void D3D::Shift3DOriginPixelPosition(float dX, float dY, float dZ) 
{
	view_shiftX += dX;
	view_shiftY += dY;
	view_shiftZ += dZ;

	g_View =
		XMMatrixLookAtRH(Eye, At, Up) *
		XMMatrixTranslation((-(float)wndWidth / 2 + view_shiftX) / pixperUnit,
		((float)wndHeight / 2 - view_shiftY) / pixperUnit,
			view_shiftZ / pixperUnit);
}

void D3D::Set3DOriginPixelPosition(float X, float Y, float Z) 
{ 
	view_shiftX = X;
	view_shiftY = Y;
	view_shiftZ = Z;

	g_View = 
		XMMatrixLookAtRH(Eye, At, Up ) * 
		XMMatrixTranslation((-(float)wndWidth/2 + view_shiftX) / pixperUnit,
							((float)wndHeight/2 - view_shiftY) / pixperUnit, 
							view_shiftZ / pixperUnit);
}

void D3D::ResetCamera(void) {

	//At vector is the origin by default
	At = XMVectorSet(0.0f, 0.0f, 0.0f, 0.0f);

	//Initialize the view matrix to default settings
	SetCamera(0.0, 0.0, CAMERADISTANCE, 0.0, 1.0, 0.0, CAMERAFOV);

	Set3DOriginPixelPosition(0, 0, 0);
}

void D3D::CleanupDevice(void) {

    if(pImmediateContext) pImmediateContext->ClearState();
	
	SafeRelease(&pd3dDevice);
	SafeRelease(&pdxgiFactory);
	SafeRelease(&pImmediateContext);

	SafeRelease(&pSwapChain);
	SafeRelease(&pRenderTargetView);

	SafeRelease(&pDepthStencilView);

	SafeRelease(&pD2DRT);
	SafeRelease(&pD2DFactory);
	SafeRelease(&pWriteFactory);
    
	Initialized = false;
	d2dInitialized = false;
}

HRESULT D3D::InitDevice(HWND hWnd) {

    HRESULT hr = S_OK;

    RECT rc;
    GetClientRect(hWnd, &rc);
    wndWidth = rc.right - rc.left;
    wndHeight = rc.bottom - rc.top;

    UINT createDeviceFlags = D3D11_CREATE_DEVICE_BGRA_SUPPORT;
#ifdef _DEBUG
    createDeviceFlags |= D3D11_CREATE_DEVICE_DEBUG;
#endif

    D3D_DRIVER_TYPE driverTypes[] = {
        D3D_DRIVER_TYPE_HARDWARE,
        D3D_DRIVER_TYPE_WARP,
        D3D_DRIVER_TYPE_REFERENCE,
    };
    UINT numDriverTypes = ARRAYSIZE(driverTypes);

    D3D_FEATURE_LEVEL featureLevels[] = {
        D3D_FEATURE_LEVEL_11_1,
        D3D_FEATURE_LEVEL_11_0,
        D3D_FEATURE_LEVEL_10_1,
        D3D_FEATURE_LEVEL_10_0,
    };
	UINT numFeatureLevels = ARRAYSIZE(featureLevels);

    for(UINT driverTypeIndex = 0; driverTypeIndex < numDriverTypes; driverTypeIndex++)
    {
        driverType = driverTypes[driverTypeIndex];
        hr = D3D11CreateDevice(nullptr, driverType, nullptr, createDeviceFlags, featureLevels, numFeatureLevels, D3D11_SDK_VERSION, &pd3dDevice, &featureLevel, &pImmediateContext);

        if(hr == E_INVALIDARG)
        {
            // DirectX 11.0 platforms will not recognize D3D_FEATURE_LEVEL_11_1 so we need to retry without it
            hr = D3D11CreateDevice(nullptr, driverType, nullptr, createDeviceFlags, &featureLevels[1], numFeatureLevels - 1, D3D11_SDK_VERSION, &pd3dDevice, &featureLevel, &pImmediateContext);
        }

        if(SUCCEEDED(hr)) break;
    }
    if(FAILED(hr)) return hr;

    // Obtain DXGI factory from device (since we used nullptr for pAdapter above)

	IDXGIDevice* dxgiDevice = nullptr;
    hr = pd3dDevice->QueryInterface(__uuidof(IDXGIDevice), reinterpret_cast<void**>(&dxgiDevice));
    if(SUCCEEDED(hr)) {

        IDXGIAdapter* adapter = nullptr;
        hr = dxgiDevice->GetAdapter(&adapter);
        if(SUCCEEDED(hr)) {
            hr = adapter->GetParent(__uuidof(IDXGIFactory1), reinterpret_cast<void**>(&pdxgiFactory));
            SafeRelease(&adapter);
        }
        SafeRelease(&dxgiDevice);
    }

    if(FAILED(hr)) return hr;

    // Create swap chain
    // DirectX 11.0 systems
	DXGI_SWAP_CHAIN_DESC sd = {};
    sd.BufferCount = 1;
    sd.BufferDesc.Width = wndWidth;
    sd.BufferDesc.Height = wndHeight;
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferDesc.RefreshRate.Numerator = 60;
    sd.BufferDesc.RefreshRate.Denominator = 1;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow = hWnd;
    sd.SampleDesc.Count = 1;
    sd.SampleDesc.Quality = 0;
    sd.Windowed = TRUE;

    hr = pdxgiFactory->CreateSwapChain(pd3dDevice, &sd, &pSwapChain);

    if(FAILED(hr)) return hr;

    // Create a render target view
    ID3D11Texture2D* pBackBuffer = nullptr;
    hr = pSwapChain->GetBuffer(0, __uuidof( ID3D11Texture2D ), reinterpret_cast<void**>(&pBackBuffer));
    if(FAILED(hr)) return hr;

    hr = pd3dDevice->CreateRenderTargetView(pBackBuffer, nullptr, &pRenderTargetView);
	SafeRelease(&pBackBuffer);
    if(FAILED(hr)) return hr;
	
    // Create depth stencil texture
	ID3D11Texture2D* pDepthStencil = nullptr;

	D3D11_TEXTURE2D_DESC descDepth = {};
    descDepth.Width = wndWidth;
    descDepth.Height = wndHeight;
    descDepth.MipLevels = 1;
    descDepth.ArraySize = 1;
    descDepth.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
    descDepth.SampleDesc.Count = 1;
    descDepth.SampleDesc.Quality = 0;
    descDepth.Usage = D3D11_USAGE_DEFAULT;
    descDepth.BindFlags = D3D11_BIND_DEPTH_STENCIL;
    descDepth.CPUAccessFlags = 0;
    descDepth.MiscFlags = 0;
    hr = pd3dDevice->CreateTexture2D(&descDepth, nullptr, &pDepthStencil);
    if(FAILED(hr)) return hr;

    // Create the depth stencil view
    D3D11_DEPTH_STENCIL_VIEW_DESC descDSV;
	ZeroMemory( &descDSV, sizeof(descDSV) );
    descDSV.Format = descDepth.Format;
    descDSV.ViewDimension = D3D11_DSV_DIMENSION_TEXTURE2D;
    descDSV.Texture2D.MipSlice = 0;
    hr = pd3dDevice->CreateDepthStencilView(pDepthStencil, &descDSV, &pDepthStencilView);
	SafeRelease(&pDepthStencil);
    if(FAILED(hr)) return hr;
    pImmediateContext->OMSetRenderTargets(1, &pRenderTargetView, pDepthStencilView);

	//Create blend state (for alpha blending, i.e. transparency)
	D3D11_BLEND_DESC blendDesc;
	ZeroMemory(&blendDesc, sizeof(D3D11_BLEND_DESC));
	blendDesc.AlphaToCoverageEnable = false;
	blendDesc.IndependentBlendEnable = false;
	blendDesc.RenderTarget[0].BlendEnable = true;
	blendDesc.RenderTarget[0].SrcBlend = D3D11_BLEND_SRC_ALPHA;
	blendDesc.RenderTarget[0].DestBlend = D3D11_BLEND_INV_SRC_ALPHA;
	blendDesc.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
	blendDesc.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ONE;
	blendDesc.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_ONE;
	blendDesc.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
	blendDesc.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;

	ID3D11BlendState * blendState;
	hr = pd3dDevice->CreateBlendState(&blendDesc, &blendState);
	if (FAILED(hr)) return hr;
	pImmediateContext->OMSetBlendState(blendState, NULL, 0xffffffff);

    // Setup the viewport
	viewPort.Width = (FLOAT)wndWidth;
	viewPort.Height = (FLOAT)wndHeight;
	viewPort.MinDepth = 0.0f;
	viewPort.MaxDepth = 1.0f;
	viewPort.TopLeftX = 0;
	viewPort.TopLeftY = 0;
    pImmediateContext->RSSetViewports(1, &viewPort);

    return S_OK;
}

HRESULT D3D::SetD2DInterop(void) 
{
	if(!Initialized) return S_FALSE;

	//Should initialize com but note it will be initialized only on the thread calling this method, so you may need to initialized it again in other routines.
	CoInitialize(nullptr);

	HRESULT hr = S_OK;

	// create the D2D factory
	D2D1_FACTORY_OPTIONS options;
	options.debugLevel = D2D1_DEBUG_LEVEL_INFORMATION;
	hr = D2D1CreateFactory(D2D1_FACTORY_TYPE_SINGLE_THREADED, options, &pD2DFactory);
	if(FAILED(hr)) return hr;

	// set up the D2D render target using the back buffer
	IDXGISurface *dxgiBackbuffer = nullptr;
	hr = pSwapChain->GetBuffer(0, IID_PPV_ARGS(&dxgiBackbuffer));
	if(FAILED(hr)) return hr;

	D2D1_RENDER_TARGET_PROPERTIES props = D2D1::RenderTargetProperties(D2D1_RENDER_TARGET_TYPE_DEFAULT, D2D1::PixelFormat(DXGI_FORMAT_UNKNOWN, D2D1_ALPHA_MODE_PREMULTIPLIED));
	hr = pD2DFactory->CreateDxgiSurfaceRenderTarget(dxgiBackbuffer, props, &pD2DRT);
	SafeRelease(&dxgiBackbuffer);
	if(FAILED(hr)) return hr;

	// create the DWrite factory
	hr = DWriteCreateFactory(DWRITE_FACTORY_TYPE_SHARED, __uuidof(pWriteFactory), (IUnknown**)(&pWriteFactory));
	if(FAILED(hr)) return hr;

	return S_OK;
}

std::string D3D::GetErrorMessage(HRESULT hr) 
{
	_com_error err(hr);
	std::wstring errMsg = err.ErrorMessage();
	std::string retString(errMsg.begin(), errMsg.end());

	return retString;
}

void D3D::DrawCBObjectBatch(std::vector<CBObjectTransform> &cbOTBatch, ObjectBufferCollection &objCol, CDO_ objColSelector) 
{
	//set object buffer and shaders to immediate context, if not already set
	objCol.SetContext(pImmediateContext, objColSelector);

	//setup constant buffer for individual transforms
    ConstantBuffer cb;
	//setup camera view and fov - these are common to all objects
	cb.mView = XMMatrixTranspose(g_View);						//the camera viewing the world xyz system - adjust using XMMatrixLookAtLH
	cb.mProjection = XMMatrixTranspose(g_Projection);			//the projection of the camera view onto the 2D screen - adjust using XMMatrixPerspectiveFovLH

	// Set primitive topology 
	pImmediateContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	for (int idx = 0; idx < cbOTBatch.size(); idx++) {

		if (cbOTBatch[idx].skip_render) continue;

		g_World = cbOTBatch[idx].Transform();
		cb.mWorld = XMMatrixTranspose(g_World);
		cb.vOutputColor = cbOTBatch[idx].Color;

		//update and draw transformed object to back buffer
		pImmediateContext->UpdateSubresource(objCol.objBuf[objColSelector]->pCB, 0, nullptr, &cb, 0, 0);
		pImmediateContext->DrawIndexed(objCol.objBuf[objColSelector]->iCount, 0, 0);
	}
}

void D3D::DrawCBObject(XMMATRIX Rotation, XMMATRIX Scaling, XMMATRIX Translation, XMFLOAT4 Color, ObjectBufferCollection &objCol, CDO_ objColSelector) 
{
	if(objColSelector < (CDO_)0) return;

	//set object buffer and shaders to immediate context, if not already set
	objCol.SetContext(pImmediateContext, objColSelector);

	//setup constant buffer for individual transforms
    ConstantBuffer cb;
	//setup camera view and fov - these are common to all objects
	cb.mView = XMMatrixTranspose(g_View);						//the camera viewing the world xyz system - adjust using XMMatrixLookAtLH
	cb.mProjection = XMMatrixTranspose(g_Projection);			//the projection of the camera view onto the 2D screen - adjust using XMMatrixPerspectiveFovLH

	g_World = Rotation*Scaling*Translation;
	cb.mWorld = XMMatrixTranspose(g_World);
	cb.vOutputColor = Color;

	// Set primitive topology 
	pImmediateContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	//update and draw transformed object to back buffer
	pImmediateContext->UpdateSubresource(objCol.objBuf[objColSelector]->pCB, 0, nullptr, &cb, 0, 0);
	pImmediateContext->DrawIndexed(objCol.objBuf[objColSelector]->iCount, 0, 0);
}

//Draw a single object using linelist primitive topology
void D3D::DrawFrameObject(XMMATRIX Rotation, XMMATRIX Scaling, XMMATRIX Translation, XMFLOAT4 Color, ObjectBufferCollection &objCol, CDO_ objColSelector)
{
	if (objColSelector < (CDO_)0) return;

	//set object buffer and shaders to immediate context, if not already set
	objCol.SetContext(pImmediateContext, objColSelector);

	//setup constant buffer for individual transforms
	ConstantBuffer cb;
	//setup camera view and fov - these are common to all objects
	cb.mView = XMMatrixTranspose(g_View);						//the camera viewing the world xyz system - adjust using XMMatrixLookAtLH
	cb.mProjection = XMMatrixTranspose(g_Projection);			//the projection of the camera view onto the 2D screen - adjust using XMMatrixPerspectiveFovLH

	g_World = Rotation * Scaling*Translation;
	cb.mWorld = XMMatrixTranspose(g_World);
	cb.vOutputColor = Color;

	// Set primitive topology 
	pImmediateContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_LINELIST);

	//update and draw transformed object to back buffer
	pImmediateContext->UpdateSubresource(objCol.objBuf[objColSelector]->pCB, 0, nullptr, &cb, 0, 0);
	pImmediateContext->DrawIndexed(objCol.objBuf[objColSelector]->iCount, 0, 0);
}

bool D3D::SaveScreenToPNG(std::string fileName, D2D1_RECT_F capture_rect) 
{
	//lock the std::mutex as this method can be called directly and interferes with rendering (e.g. call to save an image arrives just as a frame is being rendered)
	GraphicsMutex.lock();

	CoInitialize(nullptr);

	HRESULT hr = S_OK;

	//convert std::string
	std::wstring pfileName = StringtoWideString(fileName);

	ID3D11Texture2D *pDestTexture = nullptr;
	ID3D11Texture2D *pSourceTexture = nullptr;

	//Get the back buffer ... we will copy a central region of given dimensions from it
	pSwapChain->GetBuffer(0, __uuidof( ID3D11Texture2D ), reinterpret_cast<void**>(&pSourceTexture));
	
	//Make the destination 2D texture. This must have the same description as the back buffer texture, but can vary in dimensions. Set it to the required dimensions.
	D3D11_TEXTURE2D_DESC descDepth = {};
	pSourceTexture->GetDesc(&descDepth);
	descDepth.Width = (UINT)(capture_rect.right - capture_rect.left);
    descDepth.Height = (UINT)(capture_rect.bottom - capture_rect.top);
	pd3dDevice->CreateTexture2D(&descDepth, nullptr, &pDestTexture);
	
	//define the source region from which we will copy
	D3D11_BOX sourceRegion;
	sourceRegion.left = (UINT)capture_rect.left;
	sourceRegion.right = (UINT)capture_rect.right;
	sourceRegion.top = (UINT)capture_rect.top;
	sourceRegion.bottom = (UINT)capture_rect.bottom;
	sourceRegion.front = 0;
	sourceRegion.back = 1;
	
	//copy the central rectangle from source to destination
	pImmediateContext->CopySubresourceRegion( pDestTexture, 0, 0, 0, 0, pSourceTexture, 0, &sourceRegion );

	//Now save the texture containing the central part of the back buffer to the given png file -> NOTE : this needs the ScreenGrab.h file included, defined in DirectXTK_Desktop project
	hr = DirectXTK::SaveWICTextureToFile(pImmediateContext, pDestTexture, GUID_ContainerFormatPng, pfileName.c_str());

	//clean-up and return
	SafeRelease(&pDestTexture);
	SafeRelease(&pSourceTexture);

	GraphicsMutex.unlock();

	if(SUCCEEDED(hr)) return true;
	else return false;
}

HRESULT D3D::LoadScaledBitmap(std::string fileName, IWICFormatConverter **ppWICFormatConverter, D2D1_SIZE_U size) {
	
	//Return IWICFormatConverter loaded from fileName and scaled to given size.
	//If size is zero then do not rescale

	HRESULT hr = S_OK;

	//convert std::string
	std::wstring pfileName = StringtoWideString(fileName);

	//stuff needed for d2d operations

	IWICImagingFactory *pIWICFactory = nullptr;
	IWICBitmapDecoder *pDecoder = nullptr;
	IWICBitmapFrameDecode *pDecoderFrame = nullptr;
	IWICBitmapScaler *pScaler = nullptr;

	//Note : CLSID_WICImagingFactory1 for Win7 CLI !!!!! Using CLSID_WICImagingFactory will result in hr = REGDB_E_CLASSNOTREG
	//CLSID_WICImagingFactory must be used for Win32
	CoInitialize(nullptr);
	//hr = CoCreateInstance(CLSID_WICImagingFactory1, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&pIWICFactory));
	hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&pIWICFactory));
	
	//Decoder from file
	if (SUCCEEDED(hr)) hr = pIWICFactory->CreateDecoderFromFilename(pfileName.c_str(), nullptr, GENERIC_READ, WICDecodeMetadataCacheOnLoad, &pDecoder);
	
	//Get first frame from decoder
	if (SUCCEEDED(hr)) hr = pDecoder->GetFrame(0, &pDecoderFrame);
	
	//Make a scaler and scale the bitmap source (the decoder frame) or load at original size
	if (size.width && size.height) {
		if (SUCCEEDED(hr)) hr = pIWICFactory->CreateBitmapScaler(&pScaler);
		if (SUCCEEDED(hr)) hr = pScaler->Initialize(pDecoderFrame, size.width, size.height, WICBitmapInterpolationModeFant);
	}
	else {
		
		UINT bmpWidth = 0, bmpHeight = 0;

		if (SUCCEEDED(hr)) hr = pDecoderFrame->GetSize(&bmpWidth, &bmpHeight);
		if (SUCCEEDED(hr)) hr = pIWICFactory->CreateBitmapScaler(&pScaler);
		if (SUCCEEDED(hr)) hr = pScaler->Initialize(pDecoderFrame, bmpWidth, bmpHeight, WICBitmapInterpolationModeFant);
	}
		
	//Make a WIC bitmap and initialize it using the scaled decoder frame
	if (SUCCEEDED(hr)) hr = pIWICFactory->CreateFormatConverter(ppWICFormatConverter);
	if (SUCCEEDED(hr)) hr = (*ppWICFormatConverter)->Initialize(pScaler, GUID_WICPixelFormat32bppPBGRA, WICBitmapDitherTypeNone, nullptr, 0.f, WICBitmapPaletteTypeMedianCut);
	
	//clean up and return
	SafeRelease(&pScaler);
	SafeRelease(&pDecoderFrame);
	SafeRelease(&pDecoder);
	SafeRelease(&pIWICFactory);

	return hr;
}

HRESULT D3D::GetImageDimensions(std::string fileName, UINT *pbmpWidth, UINT *pbmpHeight) {

	*pbmpWidth = 0;
	*pbmpHeight = 0;

	HRESULT hr = S_OK;

	//convert std::string
	std::wstring pfileName = StringtoWideString(fileName);

	//stuff needed for d2d operations

	IWICImagingFactory *pIWICFactory = nullptr;
	IWICBitmapDecoder *pDecoder = nullptr;
	IWICBitmapFrameDecode *pDecoderFrame = nullptr;

	//Note : CLSID_WICImagingFactory1 for Win7 CLI !!!!! Using CLSID_WICImagingFactory will result in hr = REGDB_E_CLASSNOTREG
	//CLSID_WICImagingFactory must be used for Win32
	CoInitialize(nullptr);
	hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&pIWICFactory));
	
	//Decoder from file
	if (SUCCEEDED(hr)) hr = pIWICFactory->CreateDecoderFromFilename(pfileName.c_str(), nullptr, GENERIC_READ, WICDecodeMetadataCacheOnLoad, &pDecoder);
	
	//Get first frame from decoder
	if (SUCCEEDED(hr)) hr = pDecoder->GetFrame(0, &pDecoderFrame);
	
	//get dimensions in pixels of image
	if (SUCCEEDED(hr)) hr = pDecoderFrame->GetSize(pbmpWidth, pbmpHeight);

	SafeRelease(&pDecoderFrame);
	SafeRelease(&pDecoder);
	SafeRelease(&pIWICFactory);

	return hr;
}

HRESULT D3D::InitializeSinkWriter(std::string fileName, IMFSinkWriter **ppWriter, DWORD *pStreamIndex, UINT32 VIDEO_WIDTH, UINT32 VIDEO_HEIGHT, UINT32 VIDEO_FPS, int quality) {

	GUID   VIDEO_ENCODING_FORMAT = MFVideoFormat_WMV3;
	GUID   VIDEO_INPUT_FORMAT = MFVideoFormat_RGB32;
	UINT32 VIDEO_BIT_RATE;

	switch(quality) {

	case 0:
		VIDEO_BIT_RATE = VIDEO_BIT_RATE_VLOW;
		break;

	case 1:
		VIDEO_BIT_RATE = VIDEO_BIT_RATE_LOW;
		break;

	case 2:
		VIDEO_BIT_RATE = VIDEO_BIT_RATE_MEDIUM;
		break;

	case 3:
		VIDEO_BIT_RATE = VIDEO_BIT_RATE_HIGH;
		break;

	case 4:
		VIDEO_BIT_RATE = VIDEO_BIT_RATE_VHIGH;
		break;

	case 5:
		VIDEO_BIT_RATE = VIDEO_BIT_RATE_ULTRA;
		break;

	default:
		VIDEO_BIT_RATE = VIDEO_BIT_RATE_HIGH;
		break;
	}

	VIDEO_BIT_RATE = VIDEO_BIT_RATE * VIDEO_FPS / VIDEO_BIT_RATE_REFFPS;

	*ppWriter = nullptr;
    *pStreamIndex = 0;

    IMFSinkWriter   *pSinkWriter = nullptr;
    IMFMediaType    *pMediaTypeOut = nullptr;   
    IMFMediaType    *pMediaTypeIn = nullptr;   
    DWORD           streamIndex = 0;     

	//convert std::string
	std::wstring pfileName = StringtoWideString(fileName);

    HRESULT hr = MFCreateSinkWriterFromURL(pfileName.c_str(), nullptr, nullptr, &pSinkWriter);

    // Set the output media type.
    if (SUCCEEDED(hr)) hr = MFCreateMediaType(&pMediaTypeOut);
    
	if (SUCCEEDED(hr)) hr = pMediaTypeOut->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);
    
	if (SUCCEEDED(hr)) hr = pMediaTypeOut->SetGUID(MF_MT_SUBTYPE, VIDEO_ENCODING_FORMAT);
    
	if (SUCCEEDED(hr)) hr = pMediaTypeOut->SetUINT32(MF_MT_AVG_BITRATE, VIDEO_BIT_RATE);
    
	if (SUCCEEDED(hr)) hr = pMediaTypeOut->SetUINT32(MF_MT_INTERLACE_MODE, MFVideoInterlace_Progressive); 
    
	if (SUCCEEDED(hr)) hr = MFSetAttributeSize(pMediaTypeOut, MF_MT_FRAME_SIZE, VIDEO_WIDTH, VIDEO_HEIGHT);
    
	if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(pMediaTypeOut, MF_MT_FRAME_RATE, VIDEO_FPS, 1);
    
	if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(pMediaTypeOut, MF_MT_PIXEL_ASPECT_RATIO, 1, 1);
    
	if (SUCCEEDED(hr)) hr = pSinkWriter->AddStream(pMediaTypeOut, &streamIndex);

    // Set the input media type.
    if (SUCCEEDED(hr)) hr = MFCreateMediaType(&pMediaTypeIn);   
    
	if (SUCCEEDED(hr)) hr = pMediaTypeIn->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);   

    if (SUCCEEDED(hr)) hr = pMediaTypeIn->SetGUID(MF_MT_SUBTYPE, VIDEO_INPUT_FORMAT);     

    if (SUCCEEDED(hr)) hr = pMediaTypeIn->SetUINT32(MF_MT_INTERLACE_MODE, MFVideoInterlace_Progressive);   

    if (SUCCEEDED(hr)) hr = MFSetAttributeSize(pMediaTypeIn, MF_MT_FRAME_SIZE, VIDEO_WIDTH, VIDEO_HEIGHT);   

    if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(pMediaTypeIn, MF_MT_FRAME_RATE, VIDEO_FPS, 1);   

    if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(pMediaTypeIn, MF_MT_PIXEL_ASPECT_RATIO, 1, 1);   

    if (SUCCEEDED(hr)) hr = pSinkWriter->SetInputMediaType(streamIndex, pMediaTypeIn, nullptr);   

    // Tell the sink writer to start accepting data.
    if (SUCCEEDED(hr)) hr = pSinkWriter->BeginWriting();

    // Return the pointer to the caller.
    if (SUCCEEDED(hr)) {

        *ppWriter = pSinkWriter;
        (*ppWriter)->AddRef();
        *pStreamIndex = streamIndex;
    }

    SafeRelease(&pSinkWriter);
    SafeRelease(&pMediaTypeOut);
    SafeRelease(&pMediaTypeIn);
    
	return hr;
}

HRESULT D3D::WriteFrameFromFile(std::string fileName, IMFSinkWriter *pWriter, DWORD streamIndex, const LONGLONG& timeStamp, UINT32 VIDEO_WIDTH, UINT32 VIDEO_HEIGHT, UINT32 VIDEO_FPS) {

	HRESULT hr = S_OK;

	UINT64 VIDEO_FRAME_DURATION = 10 * 1000 * 1000 / VIDEO_FPS;

	IMFSample *pSample = nullptr;
    IMFMediaBuffer *pBuffer = nullptr;

	//buffer stride : pixels width * 4 bytes (since RGBA)
    const LONG cbWidth = 4 * VIDEO_WIDTH;
	//buffer size in bytes
    const DWORD cbBuffer = cbWidth * VIDEO_HEIGHT;

	//the pixel data in bytes - we will fill this shortly with data from input file
    BYTE *pData = nullptr;

	//get a wic bitmap from file, scaled to given width and height in pixels
	IWICFormatConverter *pWICBitmap = nullptr;
	D2D1_SIZE_U size = {VIDEO_WIDTH, VIDEO_HEIGHT};
	if (SUCCEEDED(hr)) hr = LoadScaledBitmap(fileName, &pWICBitmap, size);

    // Create a new memory buffer.
    if (SUCCEEDED(hr)) hr = MFCreateMemoryBuffer(cbBuffer, &pBuffer);
	
    // Lock the buffer and copy the video frame to the buffer.
    if (SUCCEEDED(hr)) hr = pBuffer->Lock(&pData, nullptr, nullptr);
	
	WICRect rclock = { 0, 0, (INT)VIDEO_WIDTH, (INT)VIDEO_HEIGHT };
	if (SUCCEEDED(hr)) hr = pWICBitmap->CopyPixels(&rclock, cbWidth, cbBuffer, pData);

    if (pBuffer) pBuffer->Unlock();
	
    // Set the data length of the buffer.
    if (SUCCEEDED(hr)) hr = pBuffer->SetCurrentLength(cbBuffer);
	
    // Create a media sample and add the buffer to the sample.
    if (SUCCEEDED(hr)) hr = MFCreateSample(&pSample);
    
	if (SUCCEEDED(hr)) hr = pSample->AddBuffer(pBuffer);
	
    // Set the time stamp and the duration.
    if (SUCCEEDED(hr)) hr = pSample->SetSampleTime(timeStamp);
	
    if (SUCCEEDED(hr)) hr = pSample->SetSampleDuration(VIDEO_FRAME_DURATION);
	
    // Send the sample to the Sink Writer. NOTE : This method fails with "Unspecified Error" if either VIDEO_WIDTH or VIDEO_HEIGHT has odd value. (!!!) Make sure to always increase these to the next nearest even value.
    if (SUCCEEDED(hr)) hr = pWriter->WriteSample(streamIndex, pSample);

    SafeRelease(&pSample);
    SafeRelease(&pBuffer);
	SafeRelease(&pWICBitmap);

    return hr;
}

//--------------------------------------------------------------------------------------------------------------------------
//
// Construct an object buffer (ConstantBuffer) with vertex and pixel shaders from input data:
//
// vertices specified as SimpleVertex which has POSITION (XMFLOAT3), COLOR (XMFLOAT4) (vertices[] with vSize = ARRAYSIZE(vertices))
// indices specifying triangles from the given vertices (indices[] with iSize = ARRAYSIZE(indices))
// vertex and index shaders are loaded from input files and must be consisten with SimpleVertex format.
// the shaders can be loaded either from a .fx (compiled output generated) or .cso file (loaded directly from compiled file)
// When loading from a .fx file the vertex shader routine must be named (the entry point) "VS" and the pixel shader (the entry point) "PS".
//
//--------------------------------------------------------------------------------------------------------------------------
HRESULT ObjectBufferCollection::NewObjectBuffer(ID3D11Device* pd3dDevice, ID3D11DeviceContext* pImmediateContext, SimpleVertex vertices[], UINT vSize, WORD indices[], UINT iSize, std::string VSFile, std::string PSFile) {
	
	if(!vSize || !iSize) return S_FALSE;

	HRESULT hr = S_OK;

	ObjectBuffer *pOB = new ObjectBuffer();

	D3D11_BUFFER_DESC bd = {};

	// Get the vertex shader blob
    ID3DBlob* pVSBlob = nullptr;
	if(GetFileTermination(VSFile).compare(".cso") == 0)							//load from a cso file
		hr = D3DReadFileToBlob(StringtoWideString(VSFile).c_str(), &pVSBlob);
	else if(GetFileTermination(VSFile).compare(".fx") == 0)						//load from a fx file - must compile it first
		hr = CompileShaderFromFile(VSFile, "VS", "vs_4_0", &pVSBlob);
	else return S_FALSE;														//input file is neither a cso nor a fx file so method call invalid
	if(FAILED(hr)) return hr;

	// Create the vertex shader
	hr = pd3dDevice->CreateVertexShader(pVSBlob->GetBufferPointer(), pVSBlob->GetBufferSize(), nullptr, &(pOB->pVertexShader));
	if(FAILED(hr)) {

		SafeRelease(&pVSBlob);
		return hr;
	}

	// Define the input layout - standard POSITION, COLOR type (hence this method is called MakeSimpleObject)
	D3D11_INPUT_ELEMENT_DESC layout[] =
    {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 }
	};
	UINT numElements = ARRAYSIZE(layout);

	// Create the input layout
	hr = pd3dDevice->CreateInputLayout(layout, numElements, pVSBlob->GetBufferPointer(), pVSBlob->GetBufferSize(), &(pOB->pVertexLayout));
	SafeRelease(&pVSBlob);

	// Set the input layout
	pImmediateContext->IASetInputLayout(pOB->pVertexLayout);

	// Compile the pixel shader
	ID3DBlob* pPSBlob = nullptr;
	if(GetFileTermination(PSFile).compare(".cso") == 0)							//load from a cso file
		hr = D3DReadFileToBlob(StringtoWideString(PSFile).c_str(), &pPSBlob);
	else if(GetFileTermination(PSFile).compare(".fx") == 0)						//load from a fx file - must compile it first
		hr = CompileShaderFromFile(PSFile, "PS", "ps_4_0", &pPSBlob);
	else hr = S_FALSE;															//input file is neither a cso nor a fx file so method call invalid

	// Create the pixel shaders
	hr = pd3dDevice->CreatePixelShader(pPSBlob->GetBufferPointer(), pPSBlob->GetBufferSize(), nullptr, &(pOB->pPixelShader));
	SafeRelease(&pPSBlob);

	// Create vertex buffer
	D3D11_SUBRESOURCE_DATA InitData = {};

	InitData.pSysMem = vertices;

    bd.Usage = D3D11_USAGE_DEFAULT;
    bd.ByteWidth = sizeof(SimpleVertex) * vSize;
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	bd.CPUAccessFlags = 0;
	hr = pd3dDevice->CreateBuffer(&bd, &InitData, &(pOB->pVertexBuffer));

	pOB->stride = sizeof(SimpleVertex);
	pOB->offset = 0;
	
	// Set vertex buffer
	pImmediateContext->IASetVertexBuffers(0, 1, &(pOB->pVertexBuffer), &pOB->stride, &pOB->offset);

	// Create index buffer - define triangle vertices in the geometric direction (anti-clockwise) when looking perpendicularly at them from the direction they should be visible.
	InitData.pSysMem = indices;

	bd.Usage = D3D11_USAGE_DEFAULT;
    bd.ByteWidth = sizeof( WORD ) * iSize;		//number of bytes contained in the indices array
    bd.BindFlags = D3D11_BIND_INDEX_BUFFER;
	bd.CPUAccessFlags = 0;
	hr =  pd3dDevice->CreateBuffer(&bd, &InitData, &(pOB->pIndexBuffer));
	
	// Set index buffer
	pImmediateContext->IASetIndexBuffer(pOB->pIndexBuffer, DXGI_FORMAT_R16_UINT, 0);

	// Set primitive topology 
	pImmediateContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	// Create the constant buffer
	bd.Usage = D3D11_USAGE_DEFAULT;
	bd.ByteWidth = sizeof(ConstantBuffer);
	bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	bd.CPUAccessFlags = 0;

    hr = pd3dDevice->CreateBuffer(&bd, nullptr, &(pOB->pCB));
	pOB->iCount = iSize;

	objBuf.push_back(pOB);

	return S_OK;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Use this method to load shaders from a .fx file and compile them to a cso file. You will need the .fx file to be in your current working directory.
//
// Once you have a cso file you can instead just load the shaders from the cso file using e.g. D3DReadFileToBlob(L"ShaderFile_VS.cso", &pVSBlob);
// For this, you will need the compiled shader files (.cso) to be in your current working directory, but not the .fx files.
//
// To call this method using e.g. CompileShaderFromFile(L"ShaderFile.fx", "VS", "vs_4_0", &pVSBlob); - this will generate a ShaderFile_VS.cso file, 
// as well as loading the shader into the Blob. The Shader entry point name must be VS in the .fx file
//
//----------------------------------------------------------------------------------------------------------------------------------------------------
HRESULT ObjectBufferCollection::CompileShaderFromFile(std::string fileName, LPCSTR pEntryPoint, LPCSTR pShaderModel, ID3DBlob** ppBlobOut) {
    
	HRESULT hr = S_OK;

    DWORD dwShaderFlags = D3DCOMPILE_ENABLE_STRICTNESS;

    ID3DBlob* pErrorBlob = nullptr;
    hr = D3DCompileFromFile( StringtoWideString(fileName).c_str(), nullptr, nullptr, pEntryPoint, pShaderModel, 
        dwShaderFlags, 0, ppBlobOut, &pErrorBlob );
    if(FAILED(hr)) {
        if(pErrorBlob)  {
            OutputDebugStringA( reinterpret_cast<const char*>(pErrorBlob->GetBufferPointer()));
            SafeRelease(&pErrorBlob);
        }
        return hr;
    }
    SafeRelease(&pErrorBlob);

    return S_OK;
}

bool D3D::MakeVideoFromFileSequence(std::string outputFile, std::vector<std::string> &inputFiles, UINT32 VIDEO_FPS, float frameScaler, int quality) {

	HRESULT hr = S_OK;

	//we will assume the input files all have the same dimensions in pixels, so we'll get the dimensions of the first file.

	if(VIDEO_FPS <= 0) return false;

	UINT64 VIDEO_FRAME_DURATION = 10 * 1000 * 1000 / VIDEO_FPS;
	UINT32 VIDEO_FRAME_COUNT = (UINT32)inputFiles.size();

	if(!VIDEO_FRAME_COUNT) return false;

	UINT32 VIDEO_WIDTH = 0;
	UINT32 VIDEO_HEIGHT = 0;

	CoInitialize(nullptr);

	hr = GetImageDimensions(inputFiles[0], &VIDEO_WIDTH, &VIDEO_HEIGHT);

	if(!VIDEO_WIDTH || !VIDEO_HEIGHT) return "Failed to make video: input files invalid.";

	if(frameScaler > 1) frameScaler = 1.0f;
	if(frameScaler < 0.25) frameScaler = 0.25f;

	VIDEO_WIDTH *=  (UINT32)frameScaler;
	VIDEO_HEIGHT *= (UINT32)frameScaler;

	//VIDEO_WIDTH and VIDEO_HEIGHT must have even values, otherwise the pWriter->WriteSample call in WriteFrameFromFile fails.
	if(VIDEO_WIDTH % 2) VIDEO_WIDTH++;
	if(VIDEO_HEIGHT % 2) VIDEO_HEIGHT++;

    if (SUCCEEDED(hr)) hr = MFStartup(MF_VERSION);
    if (SUCCEEDED(hr)) {

        IMFSinkWriter *pSinkWriter = nullptr;
        DWORD streamIndex;

        hr = InitializeSinkWriter(outputFile, &pSinkWriter, &streamIndex, VIDEO_WIDTH, VIDEO_HEIGHT, VIDEO_FPS, quality);
		
		if (SUCCEEDED(hr)) {

            // Send frames to the sink writer.
            LONGLONG timeStamp = 0;

            for (DWORD i = 0; i < VIDEO_FRAME_COUNT; ++i) {

                hr = WriteFrameFromFile(inputFiles[i], pSinkWriter, streamIndex, timeStamp, VIDEO_WIDTH, VIDEO_HEIGHT, VIDEO_FPS);
                if (FAILED(hr)) break;
                    
                timeStamp += VIDEO_FRAME_DURATION;
            }
        }

        if (SUCCEEDED(hr)) hr = pSinkWriter->Finalize();
        
		SafeRelease(&pSinkWriter);
        MFShutdown();
    }

	if (SUCCEEDED(hr)) return true;
	else return false;
}

//from mouse coordinates on screen return corresponding world position on far clipping plane (in logical units)
DBL3 D3D::Pick_FarPlane_Point(INT2 mouse)
{
	//reset world transform
	g_World = XMMatrixIdentity();

	//The world coordinates center is located at (view_shiftX - wndWidth/2, -view_shiftY + wndHeight / 2), with the y axis direction inverted.
	//Adjust mouse coordinates so that (0,0) corresponds to the world coordinates center
	mouse.x -= view_shiftX - wndWidth / 2;
	mouse.y += -view_shiftY + wndHeight / 2;

	//will obtain point on far clipping plane, hence z = 1.0 (z = 0.0 gives the near clipping plane, but this is too close to the camera to work well)
	FXMVECTOR mousePoint = { (float)mouse.x, (float)mouse.y, 1.0f, 0.0 };

	//from adjusted mouse coordinates transform to coordinates on the far clipping plane
	FXMVECTOR farPoint_fxm = XMVector3Unproject(
		mousePoint, 
		viewPort.TopLeftX, viewPort.TopLeftY, viewPort.Width, viewPort.Height, viewPort.MinDepth, viewPort.MaxDepth, 
		g_Projection, g_View, g_World);

	//farPoint now has world coordinates on the far clipping plane
	XMFLOAT4 farPoint_xm;
	XMStoreFloat4(&farPoint_xm, farPoint_fxm);

	//return far point in logical units
	return DBL3(farPoint_xm.x, farPoint_xm.y, farPoint_xm.z);
}

bool D3D::GetBitmapFromImage(std::string bitmapFile, std::vector<unsigned char>& bitmap, INT2 pixels_size)
{
	HRESULT hr = S_OK;

	IWICFormatConverter *pWICBitmap = NULL;

	//load image file by stretching/compressing it to (n.x, n.y) pixels
	D2D1_SIZE_U size = { (UINT32)pixels_size.x, (UINT32)pixels_size.y };
	if (SUCCEEDED(hr)) hr = LoadScaledBitmap(bitmapFile, &pWICBitmap, size);

	WICRect rclock = { 0, 0, pixels_size.x, pixels_size.y };

	//each pixels occupies 4 bytes : b g r a
	if (SUCCEEDED(hr)) {

		int size = pixels_size.x * pixels_size.y * 4;

		BYTE* raw_bitmap = new BYTE[size];
		hr = pWICBitmap->CopyPixels(&rclock, pixels_size.x * 4, size, raw_bitmap);

		bitmap.assign(raw_bitmap, raw_bitmap + size);

		delete[] raw_bitmap;
	}
	else {

		bitmap.clear();
	}

	SafeRelease(&pWICBitmap);

	if (SUCCEEDED(hr)) return true;
	else return false;
}
#endif