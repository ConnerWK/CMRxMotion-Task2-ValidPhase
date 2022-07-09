 #!/usr/bin/env python3
import os
import sys
import json
import argparse
import numpy as np

sys.path = [
    os.path.abspath(__file__)[:-9]
] + sys.path
# from metrics import update_dsc, update_hd95

# Relative volume difference (RVD)   1
# symmetric volume difference (SVD)  1
# volumetric overlap error (VOE)  1
# Jaccard similarity coefficient (Jaccard)  1
# Average symmetric surface distance (ASD)  1
# Root mean square symmetric surface distance (RMSD)  1
# Maximum symmetric surface distance (MSD) 1

# import numpy as np
from medpy import metric
# from metrics2 import assd, asd


def cal_rvd(pred:np.ndarray, target:np.ndarray) -> float:
    """
    calculate the relative volume difference (RVD)
    :param pred: predicted mask
    :param target: ground truth
    :return: rvd value (relative volume difference, RVD)
    """
    rvd = (np.sum(pred)/np.sum(target)-1) * 100
    return rvd


def cal_arvd(pred:np.ndarray, target:np.ndarray) -> float:
    """
    calculate the relative volume difference (RVD)
    :param pred: predicted mask
    :param target: ground truth
    :return: rvd value (relative volume difference, RVD)
    """
    rvd = np.abs((np.sum(pred)/np.sum(target)-1) * 100)
    return rvd


# def cal_asd(segmentation: np.ndarray, reference_segmentation: np.ndarray, volume_spacing: tuple = None):
#     return asd(segmentation, reference_segmentation, voxelspacing=volume_spacing, connectivity=1)


def hd95_med_version(pred:np.ndarray, true:np.ndarray, voxelspacing: tuple = (1, 1, 1)) -> float:
    """
    calculate the hd 95 distance
    :param pred:
    :param true:
    :param voxelspacing:
    :return:
    """
    hd_95 = metric.binary.hd95(pred, true, voxelspacing) # 原版的hd95,如果两者有一者为0,则报错
    return hd_95


def cal_dsc(pred: np.ndarray, gt: np.ndarray) -> float:
    """
    \计算整个3D实例的DSC系数
    :param pred: 预测实例结果
    :param gt:  金标准结果
    :return:返回两个3D实例之间的DSC系数
    """
    inter = np.sum(pred * gt)
    sum_set = np.sum(pred + gt)
    eps = np.finfo(float).eps
    dsc = 2 * inter / (sum_set + eps)
    return dsc


def update_dsc(dscs: list, img: np.ndarray, label: np.ndarray):
    """
    calculate the DSC and append it to a list with a length same with patient cases
    :param img: predicted segmentation
    :param label: the ground truth label
    :param dscs: the dsc list contains all the dscs of each patient/case
    :return: renew the dscs list
    """
    dsc = cal_dsc(img, label)
    # print(f'DSC of {calculated[-1]}: {dsc}')
    dscs.append(dsc)


def update_rvd(rvds: list, pred: np.ndarray, label: np.ndarray):
    """
    calculate the RVD and append it to a list with a length same with patient cases
    :param img: predicted segmentation
    :param label: the ground truth label
    :param rvds: the RVD list contains all the rvds of each patient/case
    :return: renew the rvds list
    """
    rvd = cal_rvd(pred, label)
    # print(f'RVD of {calculated[-1]}: {rvd}')
    rvds.append(rvd)


def update_arvd(arvds: list, pred: np.ndarray, label: np.ndarray):
    """
    calculate the aRVD and append it to a list with a length same with patient cases
    :param img: predicted segmentation
    :param label: the ground truth label
    :param arvds: the aRVD list contains all the rvds of each patient/case
    :return: renew the arvds list
    """
    arvd = cal_arvd(pred, label)
    # print(f'aRVD of {calculated[-1]}: {arvd}')
    arvds.append(arvd)


# def update_abd(abds: list, pred: np.ndarray, label: np.ndarray, volume_spacing:tuple=None):
#     """
#     calculate the ABD and append it to a list with a length same with patient cases
#     :param img: predicted segmentation
#     :param label: the ground truth label
#     :param abds: the ABD list contains all the abds of each patient/case
#     :return: renew the abds list
#     """
#     abd = cal_asd(pred, label, volume_spacing)
#     # print(f'ABD of {calculated[-1]}: {abd}')
#     abds.append(abd)


def hd95_med_version(pred:np.ndarray, true:np.ndarray, voxelspacing: tuple = (1, 1, 1)) -> float:
    """
    calculate the hd 95 distance
    :param pred:
    :param true:
    :param voxelspacing:
    :return:
    """
    hd_95 = metric.binary.hd95(pred, true, voxelspacing) # 原版的hd95,如果两者有一者为0,则报错
    return hd_95


def update_hd95(hausdorff_distance95:list, predict:np.ndarray, label:np.ndarray, volume_spacing:tuple=None):
    """
    update hausdorff_distance of each case to the hausdorff_distance list(same length with the case list)
    :param predict: the predicted segmentation mask list
    :param label: the ground truth label mtx
    :param hausdorff_distance: HD dis list will be renewed in the end of this fun
    :param volume_spacing: the original spacing of the volume data
    :return: renew hausdorff_distance
    """
    hd95 = hd95_med_version(predict, label, volume_spacing)
    hausdorff_distance95.append(hd95)


# from file_utils import nib_load, nib_load_w_spacing, unzipfile
import gzip 
import tarfile 
import zipfile 
import rarfile 
# about nii files
# import numpy as np
import nibabel as nib


'''processing gz file'''
def ungz(filename): 
    gz_file = gzip.GzipFile(filename) 
    filename = filename[:-3] # gz文件的单文件解压就是去掉 filename 后面的 .gz 
    with open(filename, "wb+") as file: 
        file.write(gz_file.read()) 
        return filename  # 这个gzip的函数需要返回值以进一步配合untar函数 

'''processing tar ball'''
def untar(filename): 
    tar = tarfile.open(filename) 
    names = tar.getnames() 
    folder_dir = '/'.join(filename.split('/')[:-1])
    # tar本身是将文件打包,解除打包会产生很多文件,因此需要建立文件夹存放 
    # if not os.path.isdir(folder_dir): 
    #     os.mkdir(folder_dir) 
    for name in names: 
        tar.extract(name, folder_dir) 
    tar.close() 

'''processing zip file'''
def unzip(filename): 
    zip_file = zipfile.ZipFile(filename) 
    folder_dir = '/'.join(filename.split('/')[:-1])
    # # 类似tar解除打包,建立文件夹存放解压的多个文件 
    # if not os.path.isdir(folder_dir): 
    #     os.mkdir(folder_dir) 
    for names in zip_file.namelist(): 
        zip_file.extract(names, folder_dir) 
    zip_file.close() 

'''processing rar file'''
def unrar(filename): 
    rar = rarfile.RarFile(filename) 
    folder_dir = '/'.join(filename.split('/')[:-1])
    # if not os.path.isdir(folder_dir): 
    #     os.mkdir(folder_dir) 
    os.chdir(folder_dir) 
    rar.extractall() 
    rar.close() 


'''unzip ziped file'''
def unzipfile(fpth):
    if '.' in fpth: 
        suffix = fpth.split('.')[-1] 
        if suffix == 'gz': 
            new_filename = ungz(fpth) 
            os.remove(fpth) 
            if new_filename.split('.')[-1] == 'tar': 
                untar(new_filename) 
                os.remove(new_filename)   
        elif suffix == 'tar': 
            untar(fpth) 
            os.remove(fpth) 
        elif suffix == 'zip': 
            unzip(fpth) 
            os.remove(fpth) 
        elif suffix == 'rar': 
            unrar(fpth) 
            os.remove(fpth) 
        return True
    else:
        return False


def nib_load(file_name):
    if not os.path.exists(file_name):
        return np.array([-1])
    proxy = nib.load(file_name)
    data = proxy.get_fdata()
    proxy.uncache()
    return data


def nib_load_w_spacing(file_name):
    if not os.path.exists(file_name):
        return np.array([-1])
    proxy = nib.load(file_name)
    data = proxy.get_fdata()
    spacing = proxy.header.get_zooms()
    proxy.uncache()
    return data, spacing


def nib_affine(file_dir):
    proxy = nib.load(file_dir)
    return proxy.affine



def get_args():
    """Set up command-line interface and get arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--submissionfile", type=str, required=True, help="Submission File")
    parser.add_argument("-g", "--goldstandard", type=str,required=True, help="Goldstandard for scoring")
    parser.add_argument("-r", "--results", type=str, required=True, default="results.json", help="Scoring results")
    return parser.parse_args()


def cal_score_tsk2(sub_folder, gt_folder, sub_fnames, gt_fnames):
    '''
        label definition
        1 (LV), 2 (MYO) and 3 (RV)
    '''
    # LV_abd, MYO_abd, RV_abd = [], [], []
    # LV_rvd, MYO_rvd, RV_rvd = [], [], []
    # LV_hd, MYO_hd, RV_hd  = [], [], [] 
    LV_dsc,  MYO_dsc,  RV_dsc  = [], [], []
    LV_hd95, MYO_hd95, RV_hd95 = [], [], []

    for sub_fname, gt_fname in zip(sub_fnames, gt_fnames):
        # paths
        sub_fpth = os.path.join(sub_folder, sub_fname)
        gt_fpth  = os.path.join(gt_folder, gt_fname)
        # load nii files
        pred_np = nib_load(sub_fpth)
        gt_np, gt_spacing = nib_load_w_spacing(gt_fpth)

        # ground truth binaryzation
        LV_gt  = np.array(gt_np==1, dtype=np.uint8)
        MYO_gt = np.array(gt_np==2, dtype=np.uint8)
        RV_gt  = np.array(gt_np==3, dtype=np.uint8)
        # prediction mask binaryzation
        LV_pred  = np.array(pred_np==1, dtype=np.uint8)
        MYO_pred = np.array(pred_np==2, dtype=np.uint8)
        RV_pred  = np.array(pred_np==3, dtype=np.uint8)
        # update metrics
        update_dsc(LV_dsc, LV_pred, LV_gt)
        update_dsc(MYO_dsc, MYO_pred, MYO_gt)
        update_dsc(RV_dsc, RV_pred, RV_gt)
        update_hd95(LV_hd95, LV_pred, LV_gt, gt_spacing)
        update_hd95(MYO_hd95, MYO_pred, MYO_gt, gt_spacing)
        update_hd95(RV_hd95, RV_pred, RV_gt, gt_spacing)

    '''calculate mean value of metrics'''
    # mean dsc values
    LV_mean_dsc = np.mean(LV_dsc) 
    MYO_mean_dsc = np.mean(MYO_dsc) 
    RV_mean_dsc = np.mean(RV_dsc) 
    # mean hd values
    LV_mean_hd95 = np.mean(LV_hd95) 
    MYO_mean_hd95 = np.mean(MYO_hd95) 
    RV_mean_hd95 = np.mean(RV_hd95) 

    scores = {
        "LV_dsc": LV_mean_dsc, 
        "MYO_dsc": MYO_mean_dsc,
        "RV_dsc": RV_mean_dsc,
        "LV_hd95": LV_mean_hd95,
        "MYO_hd95": MYO_mean_hd95,
        "RV_hd95": RV_mean_hd95
    }
    return scores


def main():
    """Main function."""
    args = get_args()

    # read files and ground truth
    subm_fname = args.submissionfile.split('/')[-1]
    gold_fname = args.goldstandard.split('/')[-1]
    sub_folder = 'submission'
    gt_folder  = 'goldstandard'

    if os.path.isdir(sub_folder):
        pass
    else:
        subm_bag_fpth = os.path.join(sub_folder, subm_fname+'.tar')
        gold_bag_fpth = os.path.join(gt_folder, gold_fname)
        os.system('mkdir submission')
        os.system('mkdir goldstandard')
        os.system(f'cp {args.submissionfile} {subm_bag_fpth}')
        os.system(f'cp {args.goldstandard} {gold_bag_fpth}')

        # unzip files
        unzipfile(subm_bag_fpth)
        unzipfile(gold_bag_fpth)
    
    sub_fnames = os.listdir(sub_folder)
    if len(sub_fnames) == 1:
        if os.path.isdir(sub_folder):
            sub_folder = os.path.join(sub_folder, sub_fnames[0])

    sub_fnames = os.listdir(sub_folder)
    gt_fnames  = os.listdir(gt_folder)
    sub_fnames.sort()
    gt_fnames.sort()
    
    # cal scores of task2
    scores = cal_score_tsk2(sub_folder, gt_folder, sub_fnames, gt_fnames)

    with open(args.results, "w") as out:
        results = {
            "submission_status": "SCORED",
            **scores
        }
        out.write(json.dumps(results))

if __name__ == "__main__":
    main()
