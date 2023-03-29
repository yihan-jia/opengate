#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import UID
import pydicom
import itk
import numpy as np
import tempfile
import datetime


def mhd_2_dicom_dose(img_mhd, dcm_output_path, ds={}, physical=True):
    """
    Parameters
    ----------
    img_mhd : itk image
        image to write in dicom dose file.
    dcm_output_path : path
        full path of the dicom file to save.
    ds : dictionary, optional
        dictionary with extra tags to write (e.g. a reference dicom object). The default is {}.
    physical : str, optional
        Bool for physical dose. The default is True. If False, dose is effective.

    Returns
    -------
    None.

    """
    # get info from image
    img = itk.GetArrayFromImage(img_mhd)
    img_size = img_mhd.GetLargestPossibleRegion().GetSize()
    img_spacing = img_mhd.GetSpacing()
    img_origin = img_mhd.GetOrigin()
    now = datetime.datetime.now()
    unique_id = pydicom.uid.generate_uid()  # create a new unique UID

    # create dcm file

    # Populate required values for file meta information
    file_meta = FileMetaDataset()
    file_meta.MediaStorageSOPClassUID = UID("1.2.840.10008.5.1.4.1.1.481.2")
    file_meta.MediaStorageSOPInstanceUID = UID("1.2.3")
    file_meta.ImplementationClassUID = UID("1.2.3.4")
    file_meta.TransferSyntaxUID = "1.2.840.10008.1.2"

    # Create the FileDataset instance (initially no data elements, but file_meta
    # supplied)
    dcm_ds = FileDataset(dcm_output_path, {}, file_meta=file_meta, preamble=b"\0" * 128)

    # set dicom default fields

    dcm_ds.SpecificCharacterSet = "ISO_IR 100"
    dcm_ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.481.2"
    dcm_ds.SOPInstanceUID = (
        unique_id  # '1.2.752.243.1.1.20180817170901595.1980.23430' ###
    )
    dcm_ds.Modality = "RTDOSE"
    dcm_ds.FrameIncrementPointer = pydicom.tag.BaseTag(
        0x3004000C
    )  # That is the tag corresponding to the "GridFrameOffsetVector". All RS dose files do it like this.
    dcm_ds.InstanceCreationDate = now.strftime("%Y%m%d")  #'20171121' ###
    dcm_ds.InstanceCreationTime = now.strftime("%H%M%S")  # '120041' ###
    dcm_ds.SamplesPerPixel = 1
    dcm_ds.PhotometricInterpretation = "MONOCHROME2"
    dcm_ds.PixelRepresentation = 0
    dcm_ds.BitsAllocated = 16
    dcm_ds.BitsStored = 16
    dcm_ds.HighBit = 15
    dcm_ds.DoseUnits = "GY"
    dcm_ds.GridFrameOffsetVector = [str(c) for c in range(9)]
    # FIXME: to be overwritten with value in RD!!!
    dcm_ds.DoseGridScaling = 2.410533e-06  # Scaling factor that when multiplied by the dose grid data found in Pixel Data (7FE0,0010) Attribute of the Image Pixel Module, yields grid doses in the dose units as specified by Dose Units (3004,0002).
    # set image information
    img = np.round(img).astype("uint16")
    dcm_ds.PixelData = img.tobytes()
    dcm_ds.Rows = img_size[1]
    dcm_ds.Columns = img_size[0]
    dcm_ds.NumberOfFrames = img_size[2]
    dcm_ds.PixelSpacing = [img_spacing[1], img_spacing[0]]
    dcm_ds.SliceThickness = img_spacing[2]
    dcm_ds.ImagePositionPatient = list(img_origin)
    dcm_ds.DoseType = "PHYSICAL" if physical else "EFFECTIVE"

    # set additional fields
    for k, v in ds.items():
        print(k)
        dcm_ds[k] = v

    # Set the transfer syntax
    dcm_ds.is_little_endian = True
    dcm_ds.is_implicit_VR = True

    dcm_ds.save_as(dcm_output_path)
    print("File saved.")


if __name__ == "__main__":
    mhd_path = "/home/fava/opengate/opengate/tests/output/output_test051_rtp/threeDdoseWater.mhd"
    img_mhd = itk.imread(mhd_path)
    mhd_2_dicom_dose(
        img_mhd,
        "/home/fava/opengate/opengate/tests/output/output_test051_rtp/my_dicom.dcm",
    )
