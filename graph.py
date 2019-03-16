import matplotlib.pyplot as plt 
file1 = open("E_avg.txt","r");
a=(file1.readlines());
b=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

x=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
for i in range(0,15):
    
    b[i]=a[i].replace("'","");
    b[i]=a[i].replace("\n","");
    #[[int(y) for y in x] for x in T1]
print('[')
for i in range(0,15):

    print(''+b[i]+',');
print(']');

y=[
1.015500,
1.012000,
1.011800,
1.011500,
1.011300,
1.011000,
1.010800,
1.010300,
1.009600,
1.009000,
1.008200,
1.007300,
1.006400,
1.005000,
1.003700,
]

plt.plot(x, y) 
  
# naming the x axis 
plt.xlabel('Epoch') 
# naming the y axis 
plt.ylabel('E_average') 
  
# giving a title to my graph 
plt.title('Multi Layer Perceptron') 
  
# function to show the plot 
plt.show()
