#include "stdafx.h"
#include "BorisGraphics.h"
#if OPERATING_SYSTEM == OS_WIN

BorisGraphics* GraphicalObject::pBG = nullptr;

BorisGraphics::BorisGraphics(HWND hWnd, std::string defFont, float defFontSize) : 
	D3D(hWnd), defFont(defFont), defFontSize(defFontSize), objCol(CDO_NONE),
	ProgramStateNames(this,
		{
			VINFO(camX), VINFO(camY), VINFO(camZ), VINFO(camDistance), VINFO(fovDeg),
			VINFO(camUx), VINFO(camUy), VINFO(camUz), VINFO(view_shiftX), VINFO(view_shiftY), VINFO(view_shiftZ),
			VINFO(original_wndWidth), VINFO(original_wndHeight)
		}, {})
{	
	//Color brush
	pVariableColorBrush = nullptr;

	//Text formats
	pNormalTextFormat = nullptr;
	pBoldTextFormat = nullptr;
	pItalicTextFormat = nullptr;
	pBoldItalicTextFormat = nullptr;

	Initialize();

	//save window width and height for saving/loading program state at the current resolution
	original_wndWidth = wndWidth;
	original_wndHeight = wndHeight;
}

BorisGraphics::~BorisGraphics() {

	Cleanup();
}

void BorisGraphics::RepairObjectState(void)
{
	//view_shiftX, view_shiftY, view_shiftZ are in pixels and might need adjusting if the original resolution the file was saved at differs
	double conversion_width = 1.0; 
	double conversion_height = 1.0; 

	if (original_wndWidth) conversion_width = (double)wndWidth / original_wndWidth;
	if (original_wndHeight) conversion_height = (double)wndHeight / original_wndHeight;

	view_shiftX *= conversion_width;
	view_shiftY *= conversion_height;

	//save window width and height for saving/loading program state at the current resolution
	original_wndWidth = wndWidth;
	original_wndHeight = wndHeight;

	SetCamera(camX, camY, camZ, camUx, camUy, camUz, fovDeg);
	Set3DOriginPixelPosition(view_shiftX, view_shiftY, view_shiftZ);
}

void BorisGraphics::Cleanup(void) {

	SafeRelease(&pVariableColorBrush);

	objCol.Release();

	SafeRelease(&pNormalTextFormat);
	SafeRelease(&pBoldTextFormat);
	SafeRelease(&pItalicTextFormat);
	SafeRelease(&pBoldItalicTextFormat);
}

void BorisGraphics::Initialize(void) {

	Setup3DGraphicsResources();

	Setup2DGraphicsResources();

	ResetCamera();
}

HRESULT BorisGraphics::Setup3DGraphicsResources() {
 
	HRESULT success = S_OK;

	success &= LoadVINFile("CDO_CONE.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	success &= LoadVINFile("CDO_ARROW.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	success &= LoadVINFile("CDO_HALFARROW.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	success &= LoadVINFile("CDO_CUBE.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	success &= LoadVINFile("CDO_HALFCUBE.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	success &= LoadVINFile("CDO_CUBEFRAME.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	success &= LoadVINFile("CDO_DODECA.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	success &= LoadVINFile("CDO_AXISX.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	success &= LoadVINFile("CDO_AXISY.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	success &= LoadVINFile("CDO_AXISZ.vin", "ShaderFile_VS.cso", "ShaderFile_PS.cso", objCol);
	
	return success;
}

void BorisGraphics::Setup2DGraphicsResources(void) {
	
	pD2DRT->CreateSolidColorBrush(D2D1::ColorF(D2D1::ColorF::White), &pVariableColorBrush);
	
	//font format : normal
	CreateTextFormat(defFont, defFontSize, DWRITE_FONT_WEIGHT_NORMAL, DWRITE_FONT_STYLE_NORMAL, &pNormalTextFormat);
	pNormalTextFormat->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
	pNormalTextFormat->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_NEAR);
	pNormalTextFormat->SetWordWrapping(DWRITE_WORD_WRAPPING_NO_WRAP);					//wrapping done manually, don't need in-built one

	//font format : bold
	CreateTextFormat(defFont, defFontSize, DWRITE_FONT_WEIGHT_BOLD, DWRITE_FONT_STYLE_NORMAL, &pBoldTextFormat);
	pBoldTextFormat->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
	pBoldTextFormat->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_NEAR);
	pBoldTextFormat->SetWordWrapping(DWRITE_WORD_WRAPPING_NO_WRAP);						//wrapping done manually, don't need in-built one

	//font format : italic
	CreateTextFormat(defFont, defFontSize, DWRITE_FONT_WEIGHT_NORMAL, DWRITE_FONT_STYLE_ITALIC, &pItalicTextFormat);
	pItalicTextFormat->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
	pItalicTextFormat->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_NEAR);
	pItalicTextFormat->SetWordWrapping(DWRITE_WORD_WRAPPING_NO_WRAP);					//wrapping done manually, don't need in-built one

	//font format : bold and italic
	CreateTextFormat(defFont, defFontSize, DWRITE_FONT_WEIGHT_BOLD, DWRITE_FONT_STYLE_ITALIC, &pBoldItalicTextFormat);
	pBoldItalicTextFormat->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
	pBoldItalicTextFormat->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_NEAR);
	pBoldItalicTextFormat->SetWordWrapping(DWRITE_WORD_WRAPPING_NO_WRAP);				//wrapping done manually, don't need in-built one

	//if this is a monospaced font then the width value will be correct - save it
	CalculateMonospacedFontPixelsWidth();
	
	//calculate height
	CalculateFontPixelsHeight();
}

//Load into objCol the geometrical object specified in the given VIN file, and using the given pixel and vertex compiled shaders (give filenames)
HRESULT BorisGraphics::LoadVINFile(std::string filename, std::string VSFile, std::string PSFile, ObjectBufferCollection &objCol) {

	HRESULT success = S_OK;

	std::vector<SimpleVertex> Vertices;
	std::vector<WORD> Indices;

	//Load vertices data from tab-separated block
	std::vector<float> VLoad;
	ReadBlockData(filename, "\t", VLoad, 8);
	for(int i = 0; i < VLoad.size()/8; i++) {

		Vertices.push_back( SimpleVertex(XMFLOAT3(VLoad[i*8 + 1], VLoad[i*8 + 2], VLoad[i*8 + 3]), 
										 XMFLOAT4(VLoad[i*8 + 4], VLoad[i*8 + 5], VLoad[i*8 + 6], VLoad[i*8 + 7])
										) 
						  );
	}

	//Load indices data from tab-separated block
	ReadBlockData(filename, "\t", Indices, 3);

	if(Vertices.size() && Indices.size())
		success &= objCol.NewObjectBuffer(pd3dDevice, pImmediateContext, &Vertices[0], (UINT)Vertices.size(), &Indices[0], (UINT)Indices.size(), VSFile, PSFile);
	else return S_FALSE;

	return success;
}

HRESULT BorisGraphics::CreatePath(D2D1_POINT_2F *ppathpoints, UINT32 pointsCount, ID2D1PathGeometry **ppPathGeometry) {

	HRESULT hr = pD2DFactory->CreatePathGeometry(&(*ppPathGeometry));

	if (SUCCEEDED(hr)) {

		ID2D1GeometrySink *pSink = nullptr;

		hr = (*ppPathGeometry)->Open(&pSink);

		if (SUCCEEDED(hr)) {

			pSink->SetFillMode(D2D1_FILL_MODE_WINDING);
			pSink->BeginFigure(ppathpoints[0], D2D1_FIGURE_BEGIN_FILLED);
			pSink->AddLines(ppathpoints, pointsCount);
			pSink->EndFigure(D2D1_FIGURE_END_CLOSED);
			hr = pSink->Close();
		}

		SafeRelease(&pSink);
	}

	return hr;
}

HRESULT BorisGraphics::CreateSolidBrush(D2D1_COLOR_F color, ID2D1SolidColorBrush **ppBrush) {

	return pD2DRT->CreateSolidColorBrush(color, ppBrush);
}

HRESULT BorisGraphics::CreateTextFormat(std::string fontName, FLOAT fontSize, DWRITE_FONT_WEIGHT fontWeight, DWRITE_FONT_STYLE fontStyle, IDWriteTextFormat **ppTextFormat) {

	IDWriteFactory *pDWriteFactory;

	HRESULT hr = DWriteCreateFactory(DWRITE_FACTORY_TYPE_SHARED, __uuidof(pDWriteFactory), reinterpret_cast<IUnknown **>(&pDWriteFactory));

	wchar_t *pfontName = StringtoWCHARPointer(fontName);

	if (SUCCEEDED(hr))
		hr = pDWriteFactory->CreateTextFormat(pfontName, nullptr, fontWeight, fontStyle,
		DWRITE_FONT_STRETCH_NORMAL, fontSize, L"", ppTextFormat);

	if (SUCCEEDED(hr)) {
		// Center the text horizontally and vertically by default: can be changed later as required.
		(*ppTextFormat)->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_CENTER);
		(*ppTextFormat)->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_CENTER);
	}

	if(pfontName) delete [] pfontName;
	SafeRelease(&pDWriteFactory);

	return hr;
}

void BorisGraphics::CalculateFontPixelsHeight(void) 
{
	HRESULT hr = S_OK;

	IDWriteTextLayout *pTextLayout = nullptr;
	DWRITE_TEXT_METRICS TextMetrics;

	if(SUCCEEDED(hr)) hr = pWriteFactory->CreateTextLayout(L"A", 1, pNormalTextFormat, (FLOAT)wndWidth, (FLOAT)wndHeight, &pTextLayout);
	
	if(SUCCEEDED(hr)) hr = pTextLayout->GetMetrics(&TextMetrics);
	
	//DWRITE_TEXT_METRICS struct has values specified in DIP (device independent pixels).
	//To convert to physical pixels, use the conversion : pixels = DIP * DPI/96, where DPI is the dots per inch setting (see https://msdn.microsoft.com/en-gb/library/windows/desktop/ff684173%28v=vs.85%29.aspx)

	if(!SUCCEEDED(hr)) { 

		TextMetrics.height = defFontSize;
	}

	//set height to physical pixels
	TextMetrics.height *= dpiY/96;

	SafeRelease(&pTextLayout);

	fontPixHeight = TextMetrics.height;
}

void BorisGraphics::CalculateMonospacedFontPixelsWidth(void)
{
	HRESULT hr = S_OK;

	IDWriteTextLayout *pTextLayout = nullptr;
	DWRITE_TEXT_METRICS TextMetrics;

	if (SUCCEEDED(hr)) hr = pWriteFactory->CreateTextLayout(L"A", 1, pNormalTextFormat, (FLOAT)wndWidth, (FLOAT)wndHeight, &pTextLayout);

	if (SUCCEEDED(hr)) hr = pTextLayout->GetMetrics(&TextMetrics);

	//DWRITE_TEXT_METRICS struct has values specified in DIP (device independent pixels).
	//To convert to physical pixels, use the conversion : pixels = DIP * DPI/96, where DPI is the dots per inch setting (see https://msdn.microsoft.com/en-gb/library/windows/desktop/ff684173%28v=vs.85%29.aspx)

	//set width to physical pixels
	TextMetrics.width = TextMetrics.width * dpiX / 96;

	SafeRelease(&pTextLayout);

	monospacedfontPixWidth = TextMetrics.width;
}

//Calculate FOV required to fit the given number of logical units along the viewing width. viewDims is the view port size in pixels.
float BorisGraphics::SetCameraFOVToFitUnits(float units, D2D1_SIZE_F viewDims)
{
	//Calculate FOV required to fit an object collection within xUnits * yUnits logical display units. 

	//calculate fov required to fit mesh along y and along x directions
	//NOTE:
	//Let u be the number of logical units required to fit along the screen width. Construct an isoscelles triangle with distance d from viewing area. Let alpha be the angle between the area normal and triangle side.
	//Then FOV = 2*alpha. alpha = atan(u/2d) so FOV = 2atan(u/2d)
	//If we need to fit uv units along the viewing rect, then uv = us * r, where r is the ratio of viewing width divided by screen width.
	//Then FOV = 2*atan( u*screen_width / 2*d*viewing_width)
	
	float fov = 2*atan(units*wndWidth/(2*camDistance*viewDims.width));

	SetCameraFOV(fov*180/XM_PI);

	return fov*180/XM_PI;
}

void BorisGraphics::RotateCameraAboutOrigin(float dAzim, float dPolar) {

	//form unit camera vector
	float cx = camX / camDistance;
	float cy = camY / camDistance;
	float cz = camZ / camDistance;

	//movement along camera up vector
	float cxnew = cx - dPolar*camUx;
	float cynew = cy - dPolar*camUy;
	float cznew = cz - dPolar*camUz;

	//need to renormalize
	NormalizeVector(cxnew, cynew, cznew);

	//form new camera up vector obtained from old camera up (rotate u by same amount and in the same plane as new camera vector is obtained by rotation of old camera vector)
	float dp1 = -(cxnew*camUx + cynew*camUy + cznew*camUz);
	float dp2 = +(cxnew*cx + cynew*cy + cznew*cz);
	camUx = cx*dp1 + camUx*dp2;
	camUy = cy*dp1 + camUy*dp2;
	camUz = cz*dp1 + camUz*dp2;

	//movement along w vector, where w = c x u
	cx = cxnew; cy = cynew; cz = cznew;

	float wx = cy*camUz - cz*camUy;
	float wy = cz*camUx - cx*camUz;
	float wz = cx*camUy - cy*camUx;

	cxnew = cx - dAzim*wx;
	cynew = cy - dAzim*wy;
	cznew = cz - dAzim*wz;

	//need to renormalize
	NormalizeVector(cxnew, cynew, cznew);

	//form new camera transverse vector obtained from old camera transverse vector (rotate w by same amount and in the same plane as new camera vector is obtained by rotation of old camera vector)
	dp1 = -(cxnew*wx + cynew*wy + cznew*wz);
	dp2 = +(cxnew*cx + cynew*cy + cznew*cz);
	wx = cx*dp1 + wx*dp2;
	wy = cy*dp1 + wy*dp2;
	wz = cz*dp1 + wz*dp2;

	//set new camera up vector, noting that if w = c x u, then u = w x c (all unit vectors)
	camUx = wy*cznew - wz*cynew;
	camUy = wz*cxnew - wx*cznew;
	camUz = wx*cynew - wy*cxnew;
	NormalizeVector(camUx, camUy, camUz);

	//Set new camera up
	SetCameraUp(camUx, camUy, camUz);

	//set new camera position
	SetCameraPosition(camDistance * cxnew, camDistance * cynew, camDistance * cznew);
}

void BorisGraphics::AdjustCameraDistanceFromOrigin(float dZ) {

	SetCameraDistance(camDistance + dZ);
}

void BorisGraphics::RotateCameraView(float dAngle) {

	//rotate camera up vector about the camera position vector
	//form unit camera vector
	float cx = camX / camDistance;
	float cy = camY / camDistance;
	float cz = camZ / camDistance;

	//w = c x u
	float wx = cy*camUz - cz*camUy;
	float wy = cz*camUx - cx*camUz;
	float wz = cx*camUy - cy*camUx;

	camUx = wx * sin(dAngle) + camUx * cos(dAngle);
	camUy = wy * sin(dAngle) + camUy * cos(dAngle);
	camUz = wz * sin(dAngle) + camUz * cos(dAngle);

	//Set new camera up
	SetCameraUp(camUx, camUy, camUz);
}

void BorisGraphics::ResetTransformation(void) {

	pD2DRT->SetTransform(D2D1::Matrix3x2F::Rotation(0, D2D1::Point2F(0,0)));
	pD2DRT->SetTransform(D2D1::Matrix3x2F::Scale(1, 1, D2D1::Point2F(0, 0)));
	pD2DRT->SetTransform(D2D1::Matrix3x2F::Translation(0, 0));
}

void BorisGraphics::FillRectangle(D2D1_RECT_F rect, D2D1_COLOR_F color) {

	BeginD2DDraw();

	pVariableColorBrush->SetColor(color);
	pD2DRT->FillRectangle(&rect, pVariableColorBrush);
	
	EndD2DDraw();
}

void BorisGraphics::FillPromptRectangle(D2D1_RECT_F rect) {

	pVariableColorBrush->SetColor(D2D1::ColorF(D2D1::ColorF::Gray));
	pD2DRT->FillRectangle(&rect, pVariableColorBrush);
}

void BorisGraphics::DrawTextLine(std::string text, D2D1_RECT_F textRect, FormatSpecifier fs)
{
	BeginD2DDraw();

	//background color
	pVariableColorBrush->SetColor(fs.bgrndColor);

	switch(fs.textOutline) {

	case TOO_SQUARE:
		//draw text background and square rectangle outline
		pD2DRT->FillRectangle(&textRect, pVariableColorBrush);
		pVariableColorBrush->SetColor(fs.textColor);
		pD2DRT->DrawRectangle(&textRect, pVariableColorBrush);
		break;

	case TOO_ROUND:
		{

			D2D1_ROUNDED_RECT roundedTextRect;
			roundedTextRect.rect = textRect;
			roundedTextRect.radiusX = 5.0f;
			roundedTextRect.radiusY = 5.0f;

			//draw text background and rounded rectangle outline
			pD2DRT->FillRoundedRectangle(&roundedTextRect, pVariableColorBrush);
			pVariableColorBrush->SetColor(fs.textColor);
			pD2DRT->DrawRoundedRectangle(&roundedTextRect, pVariableColorBrush);
		}
		break;

	default:
	
		//draw text background, no outline
		pD2DRT->FillRectangle(&textRect, pVariableColorBrush);
		pVariableColorBrush->SetColor(fs.textColor);
		break;
	}

	wchar_t *ptext = StringtoWCHARPointer(text);

	//draw the text
	if(!fs.italic && !fs.bold) {
		pNormalTextFormat->SetTextAlignment(fs.textAlignment);
		pNormalTextFormat->SetParagraphAlignment(fs.paragraphAlignment);
		pD2DRT->DrawText(ptext, (UINT32)text.length(), pNormalTextFormat, &textRect, pVariableColorBrush);
	}
	else if(fs.italic && !fs.bold) {
		pItalicTextFormat->SetTextAlignment(fs.textAlignment);
		pItalicTextFormat->SetParagraphAlignment(fs.paragraphAlignment);
		pD2DRT->DrawText(ptext, (UINT32)text.length(), pItalicTextFormat, &textRect, pVariableColorBrush);
	}
	else if(fs.bold && !fs.italic) {
		pBoldTextFormat->SetTextAlignment(fs.textAlignment);
		pBoldTextFormat->SetParagraphAlignment(fs.paragraphAlignment);
		pD2DRT->DrawText(ptext, (UINT32)text.length(), pBoldTextFormat, &textRect, pVariableColorBrush);
	}
	else {
		pBoldItalicTextFormat->SetTextAlignment(fs.textAlignment);
		pBoldItalicTextFormat->SetParagraphAlignment(fs.paragraphAlignment);
		pD2DRT->DrawText(ptext, (UINT32)text.length(), pBoldItalicTextFormat, &textRect, pVariableColorBrush);
	}
	
	if(ptext) delete [] ptext;

	EndD2DDraw();
}

void BorisGraphics::DrawSimpleTextLine(std::string text, FLT2 position, D2D1_COLOR_F textColor, IDWriteTextFormat** ppTextFormat)
{
	BeginD2DDraw();

	pVariableColorBrush->SetColor(textColor);

	wchar_t *ptext = StringtoWCHARPointer(text);

	//draw the text
	(*ppTextFormat)->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
	(*ppTextFormat)->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_NEAR);
	
	D2D1_RECT_F textRect = { position.x, position.y, position.x, position.y };
	pD2DRT->DrawText(ptext, (UINT32)text.length(), (*ppTextFormat), &textRect, pVariableColorBrush);

	if (ptext) delete[] ptext;

	EndD2DDraw();
}

//draw a batch of objects with precomputed transforms (CBObjectTransform). Draw using selected object from ObjectBufferCollection objCol
void BorisGraphics::DrawCBObjectBatch(std::vector<CBObjectTransform>& transformBatch, CDO_ objColSelector)
{
	D3D::DrawCBObjectBatch(transformBatch, objCol, objColSelector);
}

void BorisGraphics::DrawCoordinateSystem(D2D1_RECT_F spaceRect) {

	//Set default FOV
	float _fovDeg = fovDeg;
	SetCameraFOV(CAMERAFOV);

	//Set fixed camera distance
	float _camDistance = camDistance;
	SetCameraDistance(COORDSYSTCAMDIST);

	//Display the coordinate system in the bottom left-hand side corner
	float _view_shiftX = view_shiftX;
	float _view_shiftY = view_shiftY;
	float _view_shiftZ = view_shiftZ;
	Set3DOriginPixelPosition(spaceRect.left + COORDSYSTLENGTHPIXELS, spaceRect.bottom - COORDSYSTLENGTHPIXELS);

	//Y axis
	D3D::DrawCBObject(XMMatrixIdentity(), XMMatrixIdentity(), XMMatrixTranslation(0, 0, 0), XMFLOAT4(1, 0.8, 0, 1), objCol, CDO_AXISY);
	//Z axis
	D3D::DrawCBObject(XMMatrixRotationX(XM_PI/2), XMMatrixIdentity(), XMMatrixTranslation(0, 0, 0), XMFLOAT4(0, 1, 0, 1), objCol, CDO_AXISZ);
	//X Axis
	D3D::DrawCBObject(XMMatrixRotationZ(-XM_PI/2), XMMatrixIdentity(), XMMatrixTranslation(0, 0, 0), XMFLOAT4(1, 0, 0, 1), objCol, CDO_AXISX);

	//Go back to previous settings
	SetCameraFOV(_fovDeg);
	SetCameraDistance(_camDistance);
	Set3DOriginPixelPosition(_view_shiftX, _view_shiftY, _view_shiftZ);
}

float BorisGraphics::GetFontStringPixelsWidth(const std::string& str, const FormatSpecifier& fs) 
{
	HRESULT hr = S_OK;

	wchar_t *str_w = StringtoWCHARPointer(str);

	IDWriteTextLayout *pTextLayout = nullptr;
	DWRITE_TEXT_METRICS TextMetrics;

	//Note using (FLOAT)wndWidth*1000 : don't want the screen size to restrict the length of text I can measure.
	if(!fs.bold && !fs.italic) hr = pWriteFactory->CreateTextLayout(str_w, (UINT32)str.length(), pNormalTextFormat, (FLOAT)wndWidth*1000, (FLOAT)wndHeight, &pTextLayout);
	else if(fs.bold && !fs.italic) hr = pWriteFactory->CreateTextLayout(str_w, (UINT32)str.length(), pBoldTextFormat, (FLOAT)wndWidth*1000, (FLOAT)wndHeight, &pTextLayout);
	else if(fs.italic && !fs.bold) hr = pWriteFactory->CreateTextLayout(str_w, (UINT32)str.length(), pItalicTextFormat, (FLOAT)wndWidth*1000, (FLOAT)wndHeight, &pTextLayout);
	else hr = pWriteFactory->CreateTextLayout(str_w, (UINT32)str.length(), pBoldItalicTextFormat, (FLOAT)wndWidth*1000, (FLOAT)wndHeight, &pTextLayout);

	if(SUCCEEDED(hr)) hr = pTextLayout->GetMetrics(&TextMetrics);

	//DWRITE_TEXT_METRICS struct has values specified in DIP (device independent pixels).
	//To convert to physical pixels, use the conversion : pixels = DIP * DPI/96, where DPI is the dots per inch setting (see https://msdn.microsoft.com/en-gb/library/windows/desktop/ff684173%28v=vs.85%29.aspx)

	//set width to physical pixels
	TextMetrics.width = TextMetrics.width * dpiX/96;
	
	SafeRelease(&pTextLayout);

	if(str_w) delete [] str_w;

	return TextMetrics.width;
}

float BorisGraphics::GetMonospacedFontStringPixelsWidth(const std::string& str)
{
	return monospacedfontPixWidth * str.length();
}

bool BorisGraphics::SaveScreenToFile(std::string fileName, D2D1_RECT_F capture_rect)
{
	return SaveScreenToPNG(fileName, capture_rect);
}
#endif