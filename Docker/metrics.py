# Relative volume difference (RVD)   1
# symmetric volume difference (SVD)  1
# volumetric overlap error (VOE)  1
# Jaccard similarity coefficient (Jaccard)  1
# Average symmetric surface distance (ASD)  1
# Root mean square symmetric surface distance (RMSD)  1
# Maximum symmetric surface distance (MSD) 1

import numpy as np
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


# if __name__ == '__main__':
#     # input = torch.from_numpy(np.asarray([[0, 0, 0], [0, 0, 0], [0, 1, 0]]))
#     # target = torch.from_numpy(np.asarray([[0, 0, 0], [0, 0, 0], [0, 0, 1]]))
#     input = np.random.randint(0, 2, (50, 50, 50))
#     target = np.random.randint(0, 2, (50, 50, 50))
#     # input = sitk.GetImageFromArray(input)
#     # target = sitk.GetImageFromArray(target)
#     # TestHausdorffDistance(input, target)
#     #
#     # TestMeanSurfaceDistance(input, target)
#     hd = hd95_med_version(input, target, (0.625, 0.625, 1.5))
