#!/usr/bin/env python3
"""Validate prediction file.
Prediction files between Task 1 and 2 are pretty much exactly
the same, with the exception of one column name, where:
    - Task 1: done
    - Task 2: empty
"""
import os   
import sys
import json
import argparse
sys.path = [
    '/usr/local/bin/',
] + sys.path
# COLNAMES = {
#     "1": ['participant', 'was_preterm', 'probability'],
#     "2": ['participant', 'was_early_preterm', 'probability']
# }
from file_utils import *


def get_args():
    parser = argparse.ArgumentParser()
    """Set up command-line interface and get arguments."""
    parser.add_argument("-s", "--submission_file", help="Submission File")
    parser.add_argument("-g", "--goldstandard", required=True, help="Goldstandard for scoring")
    parser.add_argument("-e", "--entity_type", required=True, help="synapse entity type downloaded")
    parser.add_argument("-r", "--results", required=True, help="validation results")
    return parser.parse_args()


def main():
    """Main function."""
    args = get_args()

    if args.submission_file is None:
        prediction_file_status = "INVALID"
        invalid_reasons = ['Expected FileEntity type but found ' + args.entity_type]
    else:
        invalid_reasons = []
        prediction_file_status = "VALIDATED"

        subm_fname = args.submission_file.split('/')[-1]
        gold_fname = args.goldstandard.split('/')[-1]
        sub_folder = 'submission'
        gt_folder  = 'goldstandard'
        subm_bag_fpth = os.path.join(sub_folder, subm_fname)
        gold_bag_fpth = os.path.join(gt_folder, gold_fname)
        os.system('mkdir submission')
        os.system('mkdir goldstandard')
        os.system(f'cp {args.submission_file} {subm_bag_fpth}')
        os.system(f'cp {args.goldstandard} {gold_bag_fpth}')

        # unzip files
        unzipfile(subm_bag_fpth)
        unzipfile(gold_bag_fpth)
        # os.system(f"tar -xf {subm_tar_fpth}")
        # os.system(f"tar -xf {gold_tar_fpth}")
        # os.system(f"rm {subm_tar_fpth}")
        # os.system(f"rm {gold_tar_fpth}")
        
        sub_fnames = os.listdir(sub_folder)
        
        if len(sub_fnames) == 1:
            if os.path.isdir(sub_folder):
                sub_folder = os.path.join(sub_folder, sub_fnames[0])
            else:
                prediction_file_status = "INVALID"
                invalid_reasons = ['Expected multiple nifty files but receive single file' + args.entity_type]
                result = {'submission_errors': "\n".join(invalid_reasons),
                    'submission_status': prediction_file_status}
                with open(args.results, 'w') as o:
                    o.write(json.dumps(result))
                return

        sub_fnames = os.listdir(sub_folder)
        gt_fnames  = os.listdir(gt_folder)

        if len(sub_fnames) != len(gt_fnames):
            prediction_file_status = "INVALID"
            invalid_reasons = ['The number of prediction samples is less than the number of gold standards' + args.entity_type]
            result = {'submission_errors': "\n".join(invalid_reasons),
                'submission_status': prediction_file_status}
            with open(args.results, 'w') as o:
                o.write(json.dumps(result))
            return

        sub_fnames.sort()
        gt_fnames.sort()

        for sub_fname, gt_fname in zip(sub_fnames, gt_fnames):
            sub_fpth = os.path.join(sub_folder, sub_fname)
            gt_fpth  = os.path.join(gt_folder, gt_fname)

            sub_np = nib_load(sub_fpth)
            gt_np  = nib_load(gt_fpth)

            if sub_np.shape != gt_np.shape:
                invalid_reasons.append(f"Shape of Prediction file {sub_fname} doesn't match with origin file.")
                prediction_file_status = "INVALID"

    result = {'submission_errors': "\n".join(invalid_reasons),
              'submission_status': prediction_file_status}
    with open(args.results, 'w') as o:
        o.write(json.dumps(result))


if __name__ == "__main__":
    main()    
