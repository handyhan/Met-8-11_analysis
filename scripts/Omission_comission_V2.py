import numpy as np
import numpy.ma as ma
import pandas as pd
import math
from MSG_read import extract_data
from datetime import datetime, timedelta
pd.set_option('display.max_columns', 500)
pd.options.mode.chained_assignment = None  

def neighbours3x3(col,line):   #retruns array of neighbours given line and column
    cols_next = np.array([col-1,col,col+1])
    cols_next = np.tile(cols_next,3)
    lines_next = np.array([line-1,line,line+1])
    lines_next = np.repeat(lines_next,3)
    return cols_next,lines_next
  
    
def get_oms_coms_at_time(time, modis_df, msg_df):
    
        start_time = time - pd.Timedelta(minutes=8)
        end_time = time + pd.Timedelta(minutes=8)    
        sub_mod = modis_df.loc[start_time:end_time]
        sub_msg = msg_df.loc[start_time:end_time]
        
        if sub_mod.empty or sub_msg.empty:
            print "No Concurent pixels at time " + time                      ## pass if no concurrent pixels
        else:
            # find the 3x3 pixel window of each MSG AF pixel ------------------------------------------
            sub_msg = sub_msg.reset_index()
            for index, row in sub_msg.iterrows():
                neighbour_list = np.array(neighbours3x3(row['PIXEL_COL'],row['PIXEL_LINE']))
                if index == 0:
                    neighbour_list_full = neighbour_list
                else:
                    neighbour_list_full = np.concatenate((neighbour_list_full,neighbour_list),axis=1)           # put all into array
                    
            neighbour_list_full = pd.DataFrame(data=neighbour_list_full.T,columns = ['PIXEL_COL','PIXEL_LINE'])
            neighbour_list_full['MSG_TAG'] = neighbour_list_full['PIXEL_COL']*10000 + neighbour_list_full['PIXEL_LINE']
            neighbour_list_full= neighbour_list_full.drop_duplicates('MSG_TAG')
            columns = ['TIME', 'LATITUDE_F',' LONGITUDE_F',' FRP_F', 'FRP_UNCERTAINTY', 'PIXEL_SIZE', 'PIXEL_VZA']
            for i in columns:       
                neighbour_list_full[i] = np.nan                                                                 # make into panda and cat with original AF pixels    
            sub_msg  = pd.concat([sub_msg,neighbour_list_full], axis=0, join='outer', ignore_index=False, sort=False)
            sub_msg = sub_msg.drop_duplicates('MSG_TAG',keep='first')
            sub_msg = sub_msg.reset_index(drop=True)
            
            # identify correct fires, omitted fires and comissioned fires ---------------------------------------------------------- 
            unique = np.unique(sub_msg['MSG_TAG'])
            mod_match = sub_mod[sub_mod['MSG_TAG'].isin(unique)]                      # modis AF pixels found in MSG AF or neighbours
            unique = np.unique(mod_match['MSG_TAG'])
            msg_correct_detect = sub_msg[sub_msg['MSG_TAG'].isin(unique)]             # MSG AF and neighbours in which there was a correct modis detected AF
            msg_comission = sub_msg[~sub_msg['MSG_TAG'].isin(unique)]                 # MSG AF pixels and neighbours where no MODIS AFs were detected (error of commision) 
            mod_no_match = sub_mod[~sub_mod['MSG_TAG'].isin(unique)]                  # modis AF not found in MSG AF or neighbours (error of omission at these MSG_TAG pixels)
            mod_no_match = mod_no_match.reset_index()
            
            
            # for AF pixels only (no neightbours)
            msg_correct_pix = msg_correct_detect.dropna(subset=['TIME']).dropna(axis='columns').reset_index(drop=True)
            msg_comission_pix = msg_comission.dropna(subset=['TIME']).dropna(axis='columns').reset_index(drop=True)
            mod_omissions = mod_no_match
            
            msg_correct_pix['STATUS_MARK'] = 1                  
            msg_comission_pix['STATUS_MARK'] = 0
            sub_msg_processed  = pd.concat([msg_correct_pix,msg_comission_pix], axis=0, join='outer', ignore_index=False, sort=False)
            return sub_msg_processed,mod_omissions
        
            ### use this if need to find the neighbouring pixels of non-MSG detected modis AF (used for more detailed analysis of error of omission using MSG Quality mark )
            """
            for index, row in mod_no_match.iterrows():
                neighbour_list = np.array(neighbours3x3(row['MSG_COL'],row['MSG_LINE']))
                if index == 0:
                    neighbour_list_full = neighbour_list
                else:
                    neighbour_list_full = np.concatenate((neighbour_list_full,neighbour_list),axis=1)   
            
            mod_no_match_neighbours = neighbour_list_full
            mod_no_match_neighbours = pd.DataFrame(data=mod_no_match_neighbours.T,columns = ['PIXEL_COL','PIXEL_LINE'])
            mod_no_match_neighbours['MSG_TAG'] = mod_no_match_neighbours['PIXEL_COL']*10000 + mod_no_match_neighbours['PIXEL_LINE']   
            mod_no_match_neighbours = mod_no_match_neighbours.drop_duplicates('MSG_TAG',keep='first')
            #mod_no_match_neighbours.plot(x='PIXEL_LINE',y='PIXEL_COL',kind='scatter')
            """

    
# parameters --------------------------
BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox" 
lon_max = 20
lon_min = -15
lat_max = -5
lat_min = -10
    
#t = np.arange(datetime(2018,6,18), datetime(2018,6,19), timedelta(minutes=15)).astype(datetime)    # specify time range here
t = np.array([datetime(2018,6,18,12,15),datetime(2018,6,18,12,15)])                                 # supply one off times
for k in range(0,len(t),1):
    
    date = t[k].strftime('%Y%m%d')
    
    modis_df = pd.read_csv(BaseDir + "/txt_data/MODIS/MODIS_1806160000_1806192300", parse_dates=True,index_col='TIME') 
    msg_df = pd.read_csv(BaseDir + "/txt_data/SEVIRI/FRP_PIXELS_" +date, parse_dates=True,index_col='TIME') 
       
    ##cut region accorting to min/max lat/lon defiend above --------------------------------------
    msg_df = msg_df[(msg_df['LATITUDE'] > lat_min) & (msg_df['LATITUDE'] < lat_max) ]
    msg_df = msg_df[(msg_df['LONGITUDE'] > lon_min) & (msg_df['LONGITUDE'] < lon_max) ]    
    modis_df = modis_df[(modis_df['LATITUDE'] > lat_min) & (modis_df['LATITUDE'] < lat_max) ]
    modis_df = modis_df[(modis_df['LONGITUDE'] > lon_min) & (modis_df['LONGITUDE'] < lon_max) ]
    modis_df['MSG_TAG'] = modis_df['MSG_COL']*10000 + modis_df['MSG_LINE']                          # combine line and colum to get unique tag for seviri pixels
    msg_df['MSG_TAG'] = msg_df['PIXEL_COL']*10000 + msg_df['PIXEL_LINE']                            
    
    
    #msg_df.FRP.hist(bins=40)    #shows hist of FRP from each sat  
    #modis_df.FRP.hist(bins=40)
    
    sub_msg_processed,mod_omissions = get_oms_coms_at_time(t[k], modis_df, msg_df )     # run omission comission func at a given time   
                                                                                         # returns array with status mar [0] = comission, [1] = correct detection
                                                                                        # AND array with omission MODIS pixels (MSG_TAG here gives MSG omissions)
    
    
            #plot_single_map(sub_mod['LATITUDE'].values[:],sub_mod['LONGITUDE'].values[:],'r',("MODIS_FIRES_"+ t[k].strftime('%Y%m%d%H%M')))
            #plot_single_map(sub_msg['LATITUDE'].values[:],sub_msg['LONGITUDE'].values[:],'b',("MSG_FIRES"+t[k].strftime('%Y%m%d%H%M')))
            #bounds = [mod_no_match_neighbours['PIXEL_LINE'].min(),mod_no_match_neighbours['PIXEL_LINE'].max(),mod_no_match_neighbours['PIXEL_COL'].min(),mod_no_match_neighbours['PIXEL_COL'].max()]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    