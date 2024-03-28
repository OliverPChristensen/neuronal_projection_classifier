## Unused functions
from shapely.geometry import Polygon, Point
import pandas as pd


def transcripts_geometry(trascript_path):

    transcripts_df = pd.read_csv('./raw_data/EXP-2_B1-1_results.txt', sep='\t', header=None)
    transcripts_df.columns = ['x','y','UMI','ID','NA']

    transcripts_IDs = transcripts_df['ID'].unique()
    transcripts_IDs = transcripts_IDs[['FP ' not in item for item in transcripts_IDs]]

    # Create a list to store Shapely Point objects
    transcripts_list = {}

    # Iterate over the DataFrame and create Point objects
    for i, transcripts_ID in enumerate(transcripts_IDs):

        point_list = []
        transcripts_subset = transcripts_df[transcripts_df['ID'] == transcripts_ID]

        for index, row in transcripts_subset.iterrows():
            point = Point(row['x'], row['y'])
            point_list.append(point)

        
        print(f"trancsript {i}/100")
        

        transcripts_list[transcripts_ID] = point_list

        
    
    
    return transcripts

def create_count_matrix(rois, transcripts_path):

    transcripts_df = pd.read_csv('./raw_data/EXP-2_B1-1_results.txt', sep='\t', header=None)
    transcripts_df.columns = ['x','y','UMI','ID','NA']

    transcripts_IDs = transcripts_df['ID'].unique()
    transcripts_IDs = transcripts_IDs[['FP ' not in item for item in transcripts_IDs]]
    
    count_matrix = pd.DataFrame(0, index=rois.keys(), columns=transcripts_IDs)
    for roi in rois:
        polygon = Polygon(rois[roi])
        # Iterate over the DataFrame and create Point objects
        for i, transcripts_ID in enumerate(transcripts_IDs):

            transcripts_subset = transcripts_df[transcripts_df['ID'] == transcripts_ID]

            count = 0

            for index, row in transcripts_subset.iterrows():
                point = Point(row['x'], row['y'])

                if polygon.intersection(point):
                    count += 1
                
            
            count_matrix.loc[roi,transcripts_ID]

            
            print(f"trancsript {i}/100")


    return

