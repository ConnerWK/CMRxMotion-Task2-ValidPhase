 #!/usr/bin/env python3
import os
import sys
import json
import argparse
import numpy as np

sys.path = [
    '/usr/local/bin/'
] + sys.path
from metrics import update_dsc, update_hd95
from file_utils import nib_load, nib_load_w_spacing, unzipfile


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
        subm_bag_fpth = os.path.join(sub_folder, subm_fname)
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
