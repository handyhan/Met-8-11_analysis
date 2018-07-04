# -*- coding: utf-8 -*-
"""
Created on Thu May 24 12:32:27 2018

@author: Hannah.N
"""


def regression():
    #split into test and training
    random_split = np.random.randint(fires_m8['summed_FRP'].size, size=(fires_m8['summed_FRP'].size/2))
    
    fires_m8_train = fires_m8['summed_FRP'].values[random_split]
    fires_m8_test = fires_m8.drop(fires_m8.index[random_split])['summed_FRP']
    
    fires_m11_train = fires_m11['summed_FRP'].values[random_split]
    fires_m11_test = fires_m11.drop(fires_m11.index[random_split])['summed_FRP']
    
    regr =  linear_model.LinearRegression()
    regr.fit(fires_m8_train, fires_m11_train)
    fire_m11_pred = regr.predict(fires_m8_test)
    
    print('Coefficients: \n', regr.coef_)
    # The mean squared error
    print("Mean squared error: %.2f"
          % mean_squared_error(fires_m11_test, fire_m11_pred))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.2f' % r2_score(fires_m11_test, fire_m11_pred))
    
    # Plot outputs
    plt.scatter(fires_m8_test, fires_m11_test,  color='black')
    plt.plot(fires_m8_test, fire_m11_pred, color='blue', linewidth=3)
    
    plt.xticks(())
    plt.yticks(())
    
    plt.show()




#data = np.load('./msg_lc_v2018.npz')
#lc_type_11 = data['msg_lc_type'][0,:,:]
#x_position = np.arange(0,3712,1)
#print x_position
#x_position = np.tile(x_position,3712)
#y_position = np.repeat(x_position,3712)

#df = pd.DataFrame({'X':pd.np.tile(x_position,3712),'Y':pd.np.repeat(x_position,3712)})


#columns = []
#for i in range(0,22,1):
    #columns.append(str(data['lc_category_list'][i])) 
    #type_no = str(data['lc_category_list'][i])
    #df[type_no] = data['msg_lc_type'][i,:,:].flatten()


#df.to_csv('landcover_map')
   
#hdf = pd.HDFStore('landmap.h5')
#df = pd.read_csv('landcover_map')
#df.drop('Unnamed: 0.1',axis=1, inplace =True)
#df.drop('Unnamed: 0',axis=1, inplace =True)



#df = pd.HDFStore('landmap.h5')
#columns = list(df)
#for i in range(0,24,1):
#    hdf.put(columns[i],df[columns[i]])
    
#hdf.close()

#hdf.put('d1', df, format='table', data_columns=True)
#hdf.close()
    
#df.close()
#df = pd.HDFStore('landmap.h5')
#df = pd.read_hdf('landmap.h5',mode='r')
