import numpy as np
import pandas as pd
import nibabel as nb
from neuromaps import datasets, images, nulls, resampling, stats
from abagen import fetch_desikan_killiany
from neuromaps.datasets import available_annotations, fetch_annotation
from neuromaps.parcellate import Parcellater
from neuromaps import transforms
np.random.seed(2514)

def define_parcellation():
  # Now let's load in the dk atlas
  # Fetch desikan killiany annotation
  dk = fetch_desikan_killiany(surface = True)
  dk_info = pd.read_csv(dk['info'])
  dk_image = dk['image']
  
  # We need to map our parcellation onto the surface
  
  # Load in the dk map data
  map_left = nb.load(dk_image[0]).agg_data()
  map_right = nb.load(dk_image[1]).agg_data()
  
  # Define labels - including adding a none label
  labels = pd.concat([pd.DataFrame({"id" : [0],"label": ["???"]}), dk_info.loc[:,["id","label"]]]).reset_index(drop = True)
  
  # Parcellate the left hemi
  parc_left = images.construct_shape_gii(map_left, labels=labels.label,
                                         intent='NIFTI_INTENT_LABEL')
  parc_left.labeltable.get_labels_as_dict()
  np.unique(parc_left.agg_data())
  
  # Parcellate the right hemi
  parc_right = images.construct_shape_gii(map_right, labels=labels.label,
                                          intent='NIFTI_INTENT_LABEL')
  parc_right.labeltable.get_labels_as_dict()
  np.unique(parc_right.agg_data())
  
  # Relabel the parcellations
  parcellation = images.relabel_gifti((parc_left, parc_right))
  
  # Check the parcellation labels
  # print(parcellation[0].labeltable.get_labels_as_dict())
  # print(parcellation[1].labeltable.get_labels_as_dict())
  
  val_lh, count_lh = np.unique(parcellation[0].agg_data(), return_counts = True)
  pd.DataFrame({})
  
  val_rh, count_rh = np.unique(parcellation[1].agg_data(), return_counts = True)
  summary_df = pd.DataFrame({"val_lh": val_lh, "count_lh": count_lh, "val_rh": val_rh, "count_rh": count_rh})
  summary_df.loc[:,"diff"] = (summary_df.count_rh - summary_df.count_lh)/((summary_df.count_rh + summary_df.count_lh)/2)*100
  summary_df
  
  # Assign hemisphere to match results labels
  dk_info["hemisphere"] = ["lh" if x == "L" else "rh" if x == "R" else "" for x in dk_info["hemisphere"]]
  dk_info["hemi_label"] = [f"{dk_info.hemisphere[i]}_{dk_info.label[i]}" for i in range(dk_info.shape[0])]
  dk_info = dk_info.loc[dk_info.structure == "cortex",:]
  dk_info.loc[:,"id"] = list(range(1,dk_info.shape[0] + 1))
  
  
  return parcellation, dk_info

def spin_test(df1,df2, measure, atlas, parcellation, nperm = 1000, method = "pearsonr"):
    
  # Merge with dk info to preserve order
  atlas_data1 = atlas.merge(df1, on = "hemi_label")
  atlas_data2 = atlas.merge(df2, on = "hemi_label")
  
  # Generate nulls
  rotated = nulls.alexander_bloch(atlas_data1[measure], atlas='fsaverage', 
                                  n_perm=nperm, seed=2514, 
                                  parcellation = parcellation)
  
  # Spin test                                
  corr, pval = stats.compare_images(atlas_data1[measure], 
                                    atlas_data2[measure],
                                    metric = method, 
                                    ignore_zero = True,
                                    return_nulls = False,
                                    nulls=rotated)
  return corr, pval
