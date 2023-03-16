import os
import pydicom


def get_RS_file(dcm_dir, ss_ref_uid):
    # ss_ref_uid = self.rp_data.ReferencedStructureSetSequence[0].ReferencedSOPInstanceUID
    print(
        "going to try to find the file with structure set with UID '{}'".format(
            ss_ref_uid
        )
    )
    nskip = 0
    ndcmfail = 0
    nwrongtype = 0
    rs_data = None
    rs_path = None
    for s in os.listdir(dcm_dir):
        if s[-4:].lower() != ".dcm":
            nskip += 1
            print("no .dcm suffix: {}".format(s))
            continue
        try:
            # print(s)
            ds = pydicom.dcmread(os.path.join(dcm_dir, s))
            dcmtype = ds.SOPClassUID.name
        except:
            ndcmfail += 1
            continue
        if dcmtype == "RT Structure Set Storage" and ss_ref_uid == ds.SOPInstanceUID:
            print("found structure set for CT: {}".format(s))
            rs_data = ds
            rs_path = os.path.join(dcm_dir, s)
            break
        else:
            nwrongtype += 1
            # print("rejected structure set for CT: {}".format(s))
            # print("because it as a wrong SOP class ID: {}".format(dcmtype))
            # print("AND/OR because it has the wrong SOP Instance UID: {} != {}".format(ds.SOPInstanceUID,ss_ref_uid))
    if rs_data is None:
        raise RuntimeError(
            "could not find structure set with UID={}; skipped {} with wrong suffix, got {} with 'dcm' suffix but pydicom could not read it, got {} with wrong class UID and/or instance UID. It could well be that this is a commissioning plan without CT and structure set data.".format(
                ss_ref_uid, nskip, ndcmfail, nwrongtype
            )
        )
    return rs_data, rs_path
