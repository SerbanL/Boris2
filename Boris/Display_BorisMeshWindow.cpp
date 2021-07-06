#include "stdafx.h"
#include "Display_BorisMeshWindow.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////MESH WINDOW

BorisMeshWindow::BorisMeshWindow(D2D1_RECT_F ratios, INT2 winId) :
	WinSpace(ratios, winId),
	ProgramStateNames(this, { VINFO(physQRep), VINFO(meshPosition), VINFO(meshName), VINFO(typeName), VINFO(mouse_mesh_info_position), VINFO(mouse_mesh_info_value), VINFO(displayMeshInfo) }, {})
{
	//font format : bold
	pBG->CreateTextFormat(FONT, MESHTEXTSIZE, DWRITE_FONT_WEIGHT_BOLD, DWRITE_FONT_STYLE_NORMAL, &pMeshWinTextFormat);
	pMeshWinTextFormat->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
	pMeshWinTextFormat->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_NEAR);
	pMeshWinTextFormat->SetWordWrapping(DWRITE_WORD_WRAPPING_NO_WRAP);						//wrapping done manually, don't need in-built one
}

BorisMeshWindow::~BorisMeshWindow()
{
	SafeRelease(&pMeshWinTextFormat);
}

void BorisMeshWindow::DrawWindow(void) 
{
	//time the time it takes to draw the mesh window : this can be slow if there's too much detail, so in that case will want to reduce the amount of detail.
	DWORD start_time_ms = GetSystemTickCount();

	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//draw coordinate system
	pBG->DrawCoordinateSystem(spaceRect);
	
	//draw the physical quantity set to be displayed in the mesh
	physQRep.DrawPhysQRep();

	std::string text = " Focused Mesh : " + physQRep.get_focused_meshName() + " - " + physQRep.get_focused_typeName();
	pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Gray), &pMeshWinTextFormat);

	text = " Minimum  : " + physQRep.get_min_value_string();
	pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + 2 * MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Blue), &pMeshWinTextFormat);

	text = " Maximum  : " + physQRep.get_max_value_string();
	pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + 3 * MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Red), &pMeshWinTextFormat);

	if (displayMeshInfo) {

		text = " Position : " + mouse_mesh_info_position;
		pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + 4 * MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Green), &pMeshWinTextFormat);

		text = " Value    : " + mouse_mesh_info_value;
		pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + 5 * MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Green), &pMeshWinTextFormat);
	}

	//draw resizing frame if needed
	DrawResizingFrame();

	//too slow to be of any use : reduce amount of detail : this should speed up the drawing next time.
	if (GetSystemTickCount() - start_time_ms > MAXDRAWTIMEALLOWED_MS) ZoomNotchUp();
}

ActionOutcome BorisMeshWindow::NewMessage(AC_ aCode, INT2 mouse, std::string data) {

	ActionOutcome actionResult;

	//First implement responses shared by all derived classes
	actionResult = NewMessage_CommonResponses(aCode, mouse, data);

	//Next implement particular responses of this derived class
	switch (aCode) {

	case AC_MOUSEMIDDLEDOWN:
		if (windowEntered) SetCursor(LoadCursor(nullptr, IDC_CROSS));		//cross cursor: indicates rotation/zoom mode
		break;

	case AC_MOUSEMIDDLEUP:
		if (windowEntered) SetDefaultCursor();
		break;

	case AC_MOUSELEFTDOWN:
		if (windowEntered) {

			if (IsDoubleClick(aCode)) {

				//Double click occured

				if (physQRep.GetMouseInfo(INT2(mouse.x, mouse.y), &meshName, &typeName)) {

					actionResult.text = meshName;
					actionResult.SetCode(winId, mouse, AO_MESHFOCUS2);
				}
			}
			else {

				SetCursor(LoadCursor(nullptr, IDC_HAND));		//hand cursor: indicates rotation/zoom mode
			}
		}
		break;

	case AC_SHIFT_MOUSELEFTDOWN:
	{
		if (windowEntered) {

			//see if mouse is hovering over a mesh
			DBL3 newMeshPosition = DBL3(1, 1, 1);
			std::string unit;
			double meshValue;

			if (physQRep.GetMouseInfo(INT2(mouse.x, mouse.y), &meshName, &typeName, &newMeshPosition, &meshValue, &unit)) {

				actionResult.text = trim(ToString(newMeshPosition, "m"), ",");
				actionResult.SetCode(winId, mouse, AO_ADDCONSOLECOORDINATE);
			}
		}
	}
	break;

	case AC_MOUSELEFTUP:
		if (windowEntered) SetDefaultCursor();
		break;

	case AC_MOUSERIGHTDOWN:
		if (windowEntered) SetCursor(LoadCursor(nullptr, IDC_CROSS));		//cross cursor: indicates rotation/zoom mode
		break;

	case AC_MOUSERIGHTUP:
		if (windowEntered) SetDefaultCursor();
		break;

	case AC_MOUSEMOVE:
	{
		if (windowEntered) {

			//rotate camera view about its own axis (x mouse movement), zoom (y mouse movement)
			if (mouseRightDown) {

				pBG->AdjustCameraDistanceFromOrigin((mouse.j - dragStart.y)*CAMERAZSENS / pBG->wndHeight);
				pBG->RotateCameraView((float)(mouse.i - dragStart.x)*CAMERAROTSENS / pBG->wndWidth);
				dragStart = mouse;
				actionResult.AddCode(winId, mouse, AO_REFRESH);

				physQRep.SetAlphaBlend();
			}

			//rotate camera about origin
			else if (mouseMiddleDown) {

				pBG->RotateCameraAboutOrigin((float)(dragStart.x - mouse.i)*ROTATIONSENS / pBG->wndWidth, (float)(dragStart.y - mouse.j)*ROTATIONSENS / pBG->wndHeight);
				dragStart = mouse;
				actionResult.AddCode(winId, mouse, AO_REFRESH);

				physQRep.SetAlphaBlend();
			}

			//shift camera for mesh display view
			else if (mouseLeftDown) {

				INT2 cam_shift = (mouse - dragStart);
				pBG->Shift3DOriginPixelPosition(cam_shift.x, cam_shift.y);
				dragStart = mouse;
				actionResult.AddCode(winId, mouse, AO_REFRESH);

				physQRep.SetAlphaBlend();
			}
			else {

				//see if mouse is hovering over a mesh
				DBL3 newMeshPosition;
				std::string unit;
				double meshValue;

				if (physQRep.GetMouseInfo(INT2(mouse.x, mouse.y), &meshName, &typeName, &newMeshPosition, &meshValue, &unit)) {

					displayMeshInfo = true;
					actionResult.AddCode(winId, mouse, AO_REFRESH);

					//yes it is - is there new information?
					if (meshPosition != newMeshPosition) {

						//get position within mesh and corresponding value to be displayed
						meshPosition = newMeshPosition;
						mouse_mesh_info_position = ToString(meshPosition, "m");
						mouse_mesh_info_value = ToString(meshValue, unit) + std::string(" (") + meshName + " - " + typeName + std::string(")");
					}
				}
				else {

					if (displayMeshInfo) actionResult.AddCode(winId, mouse, AO_REFRESH);
					displayMeshInfo = false;
				}
			}
		}
		else displayMeshInfo = false;
	}
	break;

	case AC_MOUSEWHEEL:
	{
		if (windowEntered) {

			actionResult.SetCode(winId, mouse, AO_RECALCULATEMESHDISPLAY);

			int wheelDirection = ToNum(data);

			if (wheelDirection > 0) ZoomNotchUp();
			else ZoomNotchDn();
		}
	}
	break;

	case AC_DROPFILES:
	{
		if (IsMouseInWindow(mouse)) {

			actionResult.SetCode(winId, mouse, AO_FILEDROPPEDINMESH);
			//name of the file in text
			actionResult.text = data;
		}
	}
	break;

	default:
		break;
	}

	return actionResult;
}

void BorisMeshWindow::AutoSetMeshDisplaySettings(std::vector<PhysQ> physQ)
{
	physQRep.CalculateRepresentation_AutoSettings(physQ, spaceRect);
}

void BorisMeshWindow::UpdatePhysQRep(std::vector<PhysQ> physQ)
{
	physQRep.CalculateRepresentation(physQ);
}

//image_cropping specify normalized cropping within the mesh image, as left, bottom, right, top : 0, 0 point is left, bottom of screen as far as user is concerned.
bool BorisMeshWindow::SaveMeshImage(std::string fileName, DBL4 image_cropping)
{
	//here remember D2D wants 0, 0 point as the left, top points of screen, so need to invert role of bottom and top in image_cropping.
	double left = spaceRectDims.width * image_cropping.i + spaceRect.left;
	double top = spaceRectDims.height * image_cropping.j + spaceRect.top;

	double right = spaceRectDims.width * image_cropping.k + spaceRect.left;
	double bottom = spaceRectDims.height * image_cropping.l + spaceRect.top;

	D2D1_RECT_F saveRect = D2D1::RectF(left, top, right, bottom);

	return pBG->SaveScreenToFile(fileName, saveRect);
}

bool BorisMeshWindow::SaveImage(std::string fileName, std::vector<PhysQ> physQ)
{
	//save current display settings
	PhysQRepSettings start_settings = physQRep.get_current_settings();

	//calculate centered display settings
	PhysQRepSettings settings(physQ, pBG->wndWidth, spaceRect);
	
	//set maximum detail level and zoom on displayed quantity
	settings.detail_level = DETAILVALUEMIN;
	settings.m_to_l *= UNITSSCALEBACK_ZOOMED / UNITSSCALEBACK;
	physQRep.CalculateRepresentation_NewSettings(physQ, settings);

	//now draw centered physQ only
	pBG->BeginD3DDraw();

	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//draw the physical quantity set to be displayed in the mesh
	physQRep.DrawPhysQRep();

	//save rendered display to file cropping to have only the mesh display centered
	float image_cropping_i = 0.5 * (1.0 - spaceRectDims.height / spaceRectDims.width);
	if (image_cropping_i < 0) image_cropping_i = 0;

	float image_cropping_k = 0.5 * (1.0 + spaceRectDims.height / spaceRectDims.width);
	if (image_cropping_k > 1) image_cropping_k = 1;

	float image_cropping_l = 0.5 * (1.0 - spaceRectDims.width / spaceRectDims.height);
	if (image_cropping_l < 0) image_cropping_l = 0;

	float image_cropping_j = 0.5 * (1.0 + spaceRectDims.width / spaceRectDims.height);
	if (image_cropping_j > 1) image_cropping_j = 1;

	double left = spaceRectDims.width * image_cropping_i + spaceRect.left;
	double top = spaceRectDims.height * image_cropping_l + spaceRect.top;

	double right = spaceRectDims.width * image_cropping_k + spaceRect.left;
	double bottom = spaceRectDims.height * image_cropping_j + spaceRect.top;

	D2D1_RECT_F saveRect = D2D1::RectF(left, top, right, bottom);

	bool success = pBG->SaveScreenToFile(fileName, saveRect);

	//revert to previous settings
	physQRep.CalculateRepresentation_NewSettings(physQ, start_settings);

	//redraw with previous settings
	DrawWindow();

	pBG->EndD3DDraw();

	return success;
}

void BorisMeshWindow::ZoomNotchUp(void)
{
	physQRep.adjust_detail_level(DETAILNOTCHUP);
}

void BorisMeshWindow::ZoomNotchDn(void)
{
	physQRep.adjust_detail_level(DETAILNOTCHDN);
}

//rather than changing detail level using mouse wheel, set it directly using this
void BorisMeshWindow::SetDetailLevel(double detail_level)
{
	physQRep.set_detail_level(detail_level);
}

double BorisMeshWindow::GetDetailLevel(void)
{
	return physQRep.get_detail_level();
}

//Set mesh display render threshold values for faster rendering when we have many cells
void BorisMeshWindow::SetRenderThresholds(INT3 renderthresholds)
{
	physQRep.set_display_renderthresholds(renderthresholds);
}

INT3 BorisMeshWindow::GetRenderThresholds(void)
{
	return physQRep.get_display_renderthresholds();
}


#endif