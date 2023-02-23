#pragma once

#include <wrl.h>

//all this has been copied from DirectXTK project so I don't have to include the whole thing
namespace DirectXTK 
{
	using Microsoft::WRL::ComPtr;

	static bool g_WIC2 = false;

	inline bool _IsWIC2()
	{
		return g_WIC2;
	}

	inline IWICImagingFactory* _GetWIC()
	{
		static IWICImagingFactory* s_Factory = nullptr;

		if (s_Factory)
			return s_Factory;

#if(_WIN32_WINNT >= _WIN32_WINNT_WIN8) || defined(_WIN7_PLATFORM_UPDATE)
		HRESULT hr = CoCreateInstance(
			CLSID_WICImagingFactory2,
			nullptr,
			CLSCTX_INPROC_SERVER,
			__uuidof(IWICImagingFactory2),
			(LPVOID*)&s_Factory
		);

		if (SUCCEEDED(hr))
		{
			// WIC2 is available on Windows 8 and Windows 7 SP1 with KB 2670838 installed
			g_WIC2 = true;
		}
		else
		{
			hr = CoCreateInstance(
				CLSID_WICImagingFactory1,
				nullptr,
				CLSCTX_INPROC_SERVER,
				__uuidof(IWICImagingFactory),
				(LPVOID*)&s_Factory
			);

			if (FAILED(hr))
			{
				s_Factory = nullptr;
				return nullptr;
			}
		}
#else
		HRESULT hr = CoCreateInstance(
			CLSID_WICImagingFactory,
			nullptr,
			CLSCTX_INPROC_SERVER,
			__uuidof(IWICImagingFactory),
			(LPVOID*)&s_Factory
		);

		if (FAILED(hr))
		{
			s_Factory = nullptr;
			return nullptr;
		}
#endif

		return s_Factory;
	}

	inline static DXGI_FORMAT EnsureNotTypeless(DXGI_FORMAT fmt)
	{
		// Assumes UNORM or FLOAT; doesn't use UINT or SINT
		switch (fmt)
		{
		case DXGI_FORMAT_R32G32B32A32_TYPELESS: return DXGI_FORMAT_R32G32B32A32_FLOAT;
		case DXGI_FORMAT_R32G32B32_TYPELESS:    return DXGI_FORMAT_R32G32B32_FLOAT;
		case DXGI_FORMAT_R16G16B16A16_TYPELESS: return DXGI_FORMAT_R16G16B16A16_UNORM;
		case DXGI_FORMAT_R32G32_TYPELESS:       return DXGI_FORMAT_R32G32_FLOAT;
		case DXGI_FORMAT_R10G10B10A2_TYPELESS:  return DXGI_FORMAT_R10G10B10A2_UNORM;
		case DXGI_FORMAT_R8G8B8A8_TYPELESS:     return DXGI_FORMAT_R8G8B8A8_UNORM;
		case DXGI_FORMAT_R16G16_TYPELESS:       return DXGI_FORMAT_R16G16_UNORM;
		case DXGI_FORMAT_R32_TYPELESS:          return DXGI_FORMAT_R32_FLOAT;
		case DXGI_FORMAT_R8G8_TYPELESS:         return DXGI_FORMAT_R8G8_UNORM;
		case DXGI_FORMAT_R16_TYPELESS:          return DXGI_FORMAT_R16_UNORM;
		case DXGI_FORMAT_R8_TYPELESS:           return DXGI_FORMAT_R8_UNORM;
		case DXGI_FORMAT_BC1_TYPELESS:          return DXGI_FORMAT_BC1_UNORM;
		case DXGI_FORMAT_BC2_TYPELESS:          return DXGI_FORMAT_BC2_UNORM;
		case DXGI_FORMAT_BC3_TYPELESS:          return DXGI_FORMAT_BC3_UNORM;
		case DXGI_FORMAT_BC4_TYPELESS:          return DXGI_FORMAT_BC4_UNORM;
		case DXGI_FORMAT_BC5_TYPELESS:          return DXGI_FORMAT_BC5_UNORM;
		case DXGI_FORMAT_B8G8R8A8_TYPELESS:     return DXGI_FORMAT_B8G8R8A8_UNORM;
		case DXGI_FORMAT_B8G8R8X8_TYPELESS:     return DXGI_FORMAT_B8G8R8X8_UNORM;
		case DXGI_FORMAT_BC7_TYPELESS:          return DXGI_FORMAT_BC7_UNORM;
		default:                                return fmt;
		}
	}

	inline static HRESULT CaptureTexture(_In_ ID3D11DeviceContext* pContext,
		_In_ ID3D11Resource* pSource,
		_Inout_ D3D11_TEXTURE2D_DESC& desc,
		_Inout_ ComPtr<ID3D11Texture2D>& pStaging)
	{
		if (!pContext || !pSource)
			return E_INVALIDARG;

		D3D11_RESOURCE_DIMENSION resType = D3D11_RESOURCE_DIMENSION_UNKNOWN;
		pSource->GetType(&resType);

		if (resType != D3D11_RESOURCE_DIMENSION_TEXTURE2D)
			return HRESULT_FROM_WIN32(ERROR_NOT_SUPPORTED);

		ComPtr<ID3D11Texture2D> pTexture;
		HRESULT hr = pSource->QueryInterface(__uuidof(ID3D11Texture2D), reinterpret_cast<void**>(pTexture.GetAddressOf()));
		if (FAILED(hr))
			return hr;

		assert(pTexture);

		pTexture->GetDesc(&desc);

		ComPtr<ID3D11Device> d3dDevice;
		pContext->GetDevice(d3dDevice.GetAddressOf());

		if (desc.SampleDesc.Count > 1)
		{
			// MSAA content must be resolved before being copied to a staging texture
			desc.SampleDesc.Count = 1;
			desc.SampleDesc.Quality = 0;

			ComPtr<ID3D11Texture2D> pTemp;
			hr = d3dDevice->CreateTexture2D(&desc, 0, pTemp.GetAddressOf());
			if (FAILED(hr))
				return hr;

			assert(pTemp);

			DXGI_FORMAT fmt = EnsureNotTypeless(desc.Format);

			UINT support = 0;
			hr = d3dDevice->CheckFormatSupport(fmt, &support);
			if (FAILED(hr))
				return hr;

			if (!(support & D3D11_FORMAT_SUPPORT_MULTISAMPLE_RESOLVE))
				return E_FAIL;

			for (UINT item = 0; item < desc.ArraySize; ++item)
			{
				for (UINT level = 0; level < desc.MipLevels; ++level)
				{
					UINT index = D3D11CalcSubresource(level, item, desc.MipLevels);
					pContext->ResolveSubresource(pTemp.Get(), index, pSource, index, fmt);
				}
			}

			desc.BindFlags = 0;
			desc.MiscFlags &= D3D11_RESOURCE_MISC_TEXTURECUBE;
			desc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
			desc.Usage = D3D11_USAGE_STAGING;

			hr = d3dDevice->CreateTexture2D(&desc, 0, pStaging.GetAddressOf());
			if (FAILED(hr))
				return hr;

			assert(pStaging);

			pContext->CopyResource(pStaging.Get(), pTemp.Get());
		}
		else if ((desc.Usage == D3D11_USAGE_STAGING) && (desc.CPUAccessFlags & D3D11_CPU_ACCESS_READ))
		{
			// Handle case where the source is already a staging texture we can use directly
			pStaging = pTexture;
		}
		else
		{
			// Otherwise, create a staging texture from the non-MSAA source
			desc.BindFlags = 0;
			desc.MiscFlags &= D3D11_RESOURCE_MISC_TEXTURECUBE;
			desc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
			desc.Usage = D3D11_USAGE_STAGING;

			hr = d3dDevice->CreateTexture2D(&desc, 0, pStaging.GetAddressOf());
			if (FAILED(hr))
				return hr;

			assert(pStaging);

			pContext->CopyResource(pStaging.Get(), pSource);
		}

#if defined(_XBOX_ONE) && defined(_TITLE)

		if (d3dDevice->GetCreationFlags() & D3D11_CREATE_DEVICE_IMMEDIATE_CONTEXT_FAST_SEMANTICS)
		{
			ComPtr<ID3D11DeviceX> d3dDeviceX;
			hr = d3dDevice.As(&d3dDeviceX);
			if (FAILED(hr))
				return hr;

			ComPtr<ID3D11DeviceContextX> d3dContextX;
			hr = pContext->QueryInterface(__uuidof(ID3D11DeviceContextX), reinterpret_cast<void**>(d3dContextX.GetAddressOf()));
			if (FAILED(hr))
				return hr;

			UINT64 copyFence = d3dContextX->InsertFence(0);

			while (d3dDeviceX->IsFencePending(copyFence))
			{
				SwitchToThread();
			}
		}

#endif

		return S_OK;
	}

	inline HRESULT SaveWICTextureToFile(_In_ ID3D11DeviceContext* pContext,
		_In_ ID3D11Resource* pSource,
		_In_ REFGUID guidContainerFormat,
		_In_z_ LPCWSTR fileName,
		_In_opt_ const GUID* targetFormat = nullptr,
		_In_opt_ std::function<void(IPropertyBag2*)> setCustomProps = nullptr)
	{
		if (!fileName)
			return E_INVALIDARG;

		D3D11_TEXTURE2D_DESC desc = { 0 };
		ComPtr<ID3D11Texture2D> pStaging;
		HRESULT hr = CaptureTexture(pContext, pSource, desc, pStaging);
		if (FAILED(hr))
			return hr;

		// Determine source format's WIC equivalent
		WICPixelFormatGUID pfGuid;
		bool sRGB = false;
		switch (desc.Format)
		{
		case DXGI_FORMAT_R32G32B32A32_FLOAT:            pfGuid = GUID_WICPixelFormat128bppRGBAFloat; break;
		case DXGI_FORMAT_R16G16B16A16_FLOAT:            pfGuid = GUID_WICPixelFormat64bppRGBAHalf; break;
		case DXGI_FORMAT_R16G16B16A16_UNORM:            pfGuid = GUID_WICPixelFormat64bppRGBA; break;
		case DXGI_FORMAT_R10G10B10_XR_BIAS_A2_UNORM:    pfGuid = GUID_WICPixelFormat32bppRGBA1010102XR; break; // DXGI 1.1
		case DXGI_FORMAT_R10G10B10A2_UNORM:             pfGuid = GUID_WICPixelFormat32bppRGBA1010102; break;
		case DXGI_FORMAT_B5G5R5A1_UNORM:                pfGuid = GUID_WICPixelFormat16bppBGRA5551; break;
		case DXGI_FORMAT_B5G6R5_UNORM:                  pfGuid = GUID_WICPixelFormat16bppBGR565; break;
		case DXGI_FORMAT_R32_FLOAT:                     pfGuid = GUID_WICPixelFormat32bppGrayFloat; break;
		case DXGI_FORMAT_R16_FLOAT:                     pfGuid = GUID_WICPixelFormat16bppGrayHalf; break;
		case DXGI_FORMAT_R16_UNORM:                     pfGuid = GUID_WICPixelFormat16bppGray; break;
		case DXGI_FORMAT_R8_UNORM:                      pfGuid = GUID_WICPixelFormat8bppGray; break;
		case DXGI_FORMAT_A8_UNORM:                      pfGuid = GUID_WICPixelFormat8bppAlpha; break;

		case DXGI_FORMAT_R8G8B8A8_UNORM:
			pfGuid = GUID_WICPixelFormat32bppRGBA;
			break;

		case DXGI_FORMAT_R8G8B8A8_UNORM_SRGB:
			pfGuid = GUID_WICPixelFormat32bppRGBA;
			sRGB = true;
			break;

		case DXGI_FORMAT_B8G8R8A8_UNORM: // DXGI 1.1
			pfGuid = GUID_WICPixelFormat32bppBGRA;
			break;

		case DXGI_FORMAT_B8G8R8A8_UNORM_SRGB: // DXGI 1.1
			pfGuid = GUID_WICPixelFormat32bppBGRA;
			sRGB = true;
			break;

		case DXGI_FORMAT_B8G8R8X8_UNORM: // DXGI 1.1
			pfGuid = GUID_WICPixelFormat32bppBGR;
			break;

		case DXGI_FORMAT_B8G8R8X8_UNORM_SRGB: // DXGI 1.1
			pfGuid = GUID_WICPixelFormat32bppBGR;
			sRGB = true;
			break;

		default:
			return HRESULT_FROM_WIN32(ERROR_NOT_SUPPORTED);
		}

		IWICImagingFactory* pWIC = _GetWIC();
		if (!pWIC)
			return E_NOINTERFACE;

		ComPtr<IWICStream> stream;
		hr = pWIC->CreateStream(stream.GetAddressOf());
		if (FAILED(hr))
			return hr;

		hr = stream->InitializeFromFilename(fileName, GENERIC_WRITE);
		if (FAILED(hr))
			return hr;

		ComPtr<IWICBitmapEncoder> encoder;
		hr = pWIC->CreateEncoder(guidContainerFormat, 0, encoder.GetAddressOf());
		if (FAILED(hr))
			return hr;

		hr = encoder->Initialize(stream.Get(), WICBitmapEncoderNoCache);
		if (FAILED(hr))
			return hr;

		ComPtr<IWICBitmapFrameEncode> frame;
		ComPtr<IPropertyBag2> props;
		hr = encoder->CreateNewFrame(frame.GetAddressOf(), props.GetAddressOf());
		if (FAILED(hr))
			return hr;

		if (targetFormat && memcmp(&guidContainerFormat, &GUID_ContainerFormatBmp, sizeof(WICPixelFormatGUID)) == 0 && _IsWIC2())
		{
			// Opt-in to the WIC2 support for writing 32-bit Windows BMP files with an alpha channel
			PROPBAG2 option = { 0 };
			option.pstrName = L"EnableV5Header32bppBGRA";

			VARIANT varValue;
			varValue.vt = VT_BOOL;
			varValue.boolVal = VARIANT_TRUE;
			(void)props->Write(1, &option, &varValue);
		}

		if (setCustomProps)
		{
			setCustomProps(props.Get());
		}

		hr = frame->Initialize(props.Get());
		if (FAILED(hr))
			return hr;

		hr = frame->SetSize(desc.Width, desc.Height);
		if (FAILED(hr))
			return hr;

		hr = frame->SetResolution(72, 72);
		if (FAILED(hr))
			return hr;

		// Pick a target format
		WICPixelFormatGUID targetGuid;
		if (targetFormat)
		{
			targetGuid = *targetFormat;
		}
		else
		{
			// Screenshots don’t typically include the alpha channel of the render target
			switch (desc.Format)
			{
#if (_WIN32_WINNT >= _WIN32_WINNT_WIN8) || defined(_WIN7_PLATFORM_UPDATE)
			case DXGI_FORMAT_R32G32B32A32_FLOAT:
			case DXGI_FORMAT_R16G16B16A16_FLOAT:
				if (_IsWIC2())
				{
					targetGuid = GUID_WICPixelFormat96bppRGBFloat;
				}
				else
				{
					targetGuid = GUID_WICPixelFormat24bppBGR;
				}
				break;
#endif

			case DXGI_FORMAT_R16G16B16A16_UNORM: targetGuid = GUID_WICPixelFormat48bppBGR; break;
			case DXGI_FORMAT_B5G5R5A1_UNORM:     targetGuid = GUID_WICPixelFormat16bppBGR555; break;
			case DXGI_FORMAT_B5G6R5_UNORM:       targetGuid = GUID_WICPixelFormat16bppBGR565; break;

			case DXGI_FORMAT_R32_FLOAT:
			case DXGI_FORMAT_R16_FLOAT:
			case DXGI_FORMAT_R16_UNORM:
			case DXGI_FORMAT_R8_UNORM:
			case DXGI_FORMAT_A8_UNORM:
				targetGuid = GUID_WICPixelFormat8bppGray;
				break;

			default:
				targetGuid = GUID_WICPixelFormat24bppBGR;
				break;
			}
		}

		hr = frame->SetPixelFormat(&targetGuid);
		if (FAILED(hr))
			return hr;

		if (targetFormat && memcmp(targetFormat, &targetGuid, sizeof(WICPixelFormatGUID)) != 0)
		{
			// Requested output pixel format is not supported by the WIC codec
			return E_FAIL;
		}

		// Encode WIC metadata
		ComPtr<IWICMetadataQueryWriter> metawriter;
		if (SUCCEEDED(frame->GetMetadataQueryWriter(metawriter.GetAddressOf())))
		{
			PROPVARIANT value;
			PropVariantInit(&value);

			value.vt = VT_LPSTR;
			value.pszVal = "DirectXTK";

			if (memcmp(&guidContainerFormat, &GUID_ContainerFormatPng, sizeof(GUID)) == 0)
			{
				// Set Software name
				(void)metawriter->SetMetadataByName(L"/tEXt/{str=Software}", &value);

				// Set sRGB chunk
				if (sRGB)
				{
					value.vt = VT_UI1;
					value.bVal = 0;
					(void)metawriter->SetMetadataByName(L"/sRGB/RenderingIntent", &value);
				}
			}
#if defined(_XBOX_ONE) && defined(_TITLE)
			else if (memcmp(&guidContainerFormat, &GUID_ContainerFormatJpeg, sizeof(GUID)) == 0)
			{
				// Set Software name
				(void)metawriter->SetMetadataByName(L"/app1/ifd/{ushort=305}", &value);

				if (sRGB)
				{
					// Set EXIF Colorspace of sRGB
					value.vt = VT_UI2;
					value.uiVal = 1;
					(void)metawriter->SetMetadataByName(L"/app1/ifd/exif/{ushort=40961}", &value);
				}
			}
			else if (memcmp(&guidContainerFormat, &GUID_ContainerFormatTiff, sizeof(GUID)) == 0)
			{
				// Set Software name
				(void)metawriter->SetMetadataByName(L"/ifd/{ushort=305}", &value);

				if (sRGB)
				{
					// Set EXIF Colorspace of sRGB
					value.vt = VT_UI2;
					value.uiVal = 1;
					(void)metawriter->SetMetadataByName(L"/ifd/exif/{ushort=40961}", &value);
				}
			}
#else
			else
			{
				// Set Software name
				(void)metawriter->SetMetadataByName(L"System.ApplicationName", &value);

				if (sRGB)
				{
					// Set EXIF Colorspace of sRGB
					value.vt = VT_UI2;
					value.uiVal = 1;
					(void)metawriter->SetMetadataByName(L"System.Image.ColorSpace", &value);
				}
			}
#endif
		}

		D3D11_MAPPED_SUBRESOURCE mapped;
		hr = pContext->Map(pStaging.Get(), 0, D3D11_MAP_READ, 0, &mapped);
		if (FAILED(hr))
			return hr;

		if (memcmp(&targetGuid, &pfGuid, sizeof(WICPixelFormatGUID)) != 0)
		{
			// Conversion required to write
			ComPtr<IWICBitmap> source;
			hr = pWIC->CreateBitmapFromMemory(desc.Width, desc.Height, pfGuid,
				mapped.RowPitch, mapped.RowPitch * desc.Height,
				reinterpret_cast<BYTE*>(mapped.pData), source.GetAddressOf());
			if (FAILED(hr))
			{
				pContext->Unmap(pStaging.Get(), 0);
				return hr;
			}

			ComPtr<IWICFormatConverter> FC;
			hr = pWIC->CreateFormatConverter(FC.GetAddressOf());
			if (FAILED(hr))
			{
				pContext->Unmap(pStaging.Get(), 0);
				return hr;
			}

			BOOL canConvert = FALSE;
			hr = FC->CanConvert(pfGuid, targetGuid, &canConvert);
			if (FAILED(hr) || !canConvert)
			{
				return E_UNEXPECTED;
			}

			hr = FC->Initialize(source.Get(), targetGuid, WICBitmapDitherTypeNone, 0, 0, WICBitmapPaletteTypeCustom);
			if (FAILED(hr))
			{
				pContext->Unmap(pStaging.Get(), 0);
				return hr;
			}

			WICRect rect = { 0, 0, static_cast<INT>(desc.Width), static_cast<INT>(desc.Height) };
			hr = frame->WriteSource(FC.Get(), &rect);
			if (FAILED(hr))
			{
				pContext->Unmap(pStaging.Get(), 0);
				return hr;
			}
		}
		else
		{
			// No conversion required
			hr = frame->WritePixels(desc.Height, mapped.RowPitch, mapped.RowPitch * desc.Height, reinterpret_cast<BYTE*>(mapped.pData));
			if (FAILED(hr))
				return hr;
		}

		pContext->Unmap(pStaging.Get(), 0);

		hr = frame->Commit();
		if (FAILED(hr))
			return hr;

		hr = encoder->Commit();
		if (FAILED(hr))
			return hr;

		return S_OK;
	}

}
