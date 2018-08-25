import numpy as np
import pandas as pd
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
            print "No Concurent pixels at time " + str(time)                      ## pass if no concurrent pixels
            sub_msg_processed = pd.DataFrame()
            mod_omissions = pd.DataFrame()
            mod_match= pd.DataFrame()
            return sub_msg_processed,mod_omissions,mod_match
        else:
        
            sub_mod['MSG_TAG'] = sub_mod['MSG_COL']*10000 + sub_mod['MSG_LINE']                          # combine line and colum to get unique tag for seviri pixels
            sub_msg['MSG_TAG'] = sub_msg['PIXEL_COL']*10000 + sub_msg['PIXEL_LINE']           
            
            #remove APROXIMATELY  MSG pixels not in Modis granuale
            ulat = sub_mod['LATITUDE'].max()            
            llat = sub_mod['LATITUDE'].min()
            ulon = sub_mod['LONGITUDE'].max()
            llon = sub_mod['LONGITUDE'].min()
                      
            sub_msg = sub_msg[(sub_msg['LATITUDE'] > llat) & (sub_msg['LATITUDE'] < ulat) ]                  
            sub_msg = sub_msg[(sub_msg['LONGITUDE'] > llon) & (sub_msg['LONGITUDE'] < ulon) ]
            
            """ #### if known that there is a MYD and MOD overpas, split the subsetting of latlong        
            if len(sub_mod['SAT'].unique()) > 1 :
                print sub_mod['SAT']
                print " MYD and MOD at same time"
            else:
                print  sub_mod.sort_values('LATITUDE')
            """
            
            # find the 3x3 pixel window of each MSG AF pixel ------------------------------------------ 
            
            if sub_msg.empty:                                                              # if no MSG pixels at that time in modis region pass
                sub_msg_processed = pd.DataFrame()
                mod_omissions = pd.DataFrame()
                mod_match= pd.DataFrame()
                return sub_msg_processed,mod_omissions,mod_match
            else:
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
                sub_msg = sub_msg.append(neighbour_list_full)                           # make into panda and add with original AF pixels
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
                msg_correct_pix = msg_correct_detect.dropna().reset_index(drop=True) 
                msg_comission_pix = msg_comission.dropna().reset_index(drop=True)            # can drop non AF pixe neighbours now
                mod_omissions = mod_no_match
                msg_correct_pix['STATUS_MARK'] = 1                                          # correct detection is 1 comission is 0
                msg_comission_pix['STATUS_MARK'] = 0
                sub_msg_processed  = pd.concat([msg_correct_pix,msg_comission_pix], axis=0, join='outer', ignore_index=False, sort=False)
                return sub_msg_processed, mod_omissions, mod_match
        
        
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


"""
correct_coms = pd.read_csv(BaseDir + "/TEST_DATA/correct_n_commission_full", parse_dates=True,index_col='TIME') 
omissions = pd.read_csv(BaseDir + "/TEST_DATA/omission_full", parse_dates=True,index_col='TIME') 

correct = correct_coms[correct_coms['STATUS_MARK'] == 1 ]
coms = correct_coms[correct_coms['STATUS_MARK'] == 0 ]


correct_frac = len(correct)/(float(len(correct))+len(omissions)+len(coms))
oms_frac = len(omissions)/(float(len(correct))+len(omissions)+len(coms))
coms_frac = len(coms)/(float(len(correct))+len(omissions)+len(coms))

print correct_frac,oms_frac,coms_frac

"""

#plot_single_map(sub_mod['LATITUDE'].values[:],sub_mod['LONGITUDE'].values[:],'r',("MODIS_FIRES_"+ t[k].strftime('%Y%m%d%H%M')))
#plot_single_map(sub_msg['LATITUDE'].values[:],sub_msg['LONGITUDE'].values[:],'b',("MSG_FIRES"+t[k].strftime('%Y%m%d%H%M')))
#bounds = [mod_no_match_neighbours['PIXEL_LINE'].min(),mod_no_match_neighbours['PIXEL_LINE'].max(),mod_no_match_neighbours['PIXEL_COL'].min(),mod_no_match_neighbours['PIXEL_COL'].max()]


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    