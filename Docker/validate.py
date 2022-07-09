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
    os.path.abspath(__file__)[:-12]
] + sys.path
# COLNAMES = {
#     "1": ['participant', 'was_preterm', 'probability'],
#     "2": ['participant', 'was_early_preterm', 'probability']
# }
# from file_utils import *

# about file decomposing
import gzip 
import tarfile 
import zipfile 
import rarfile 
# about nii files
import numpy as np
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
        subm_bag_fpth = os.path.join(sub_folder, subm_fname+".tar")
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
