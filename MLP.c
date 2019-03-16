#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#define N_ATTR 10  //Number of Attributes as mentioned in dataset characteristics
#define N_CLASSES 7 // Number of classes as mentioned in dataset characteristics
#define N_SAMPLES 1000 // No of training Examples
#define OP1 1	// Output value for class1
#define OP2 -1	//Output value for class2
#define DATA_FILE "Training_EXOR.txt" // training dataset file
#define LEARNING_RATE 0.2 //value of learning parameter
#define THRESHOLD 0
#define NO_ITERATIONS 15
#define TEST_SAMPLES 1000
#define LAMBDA 0.7


//long long long long long long long long long long long long float Calculate_Weighted_Sum(float input[],int size,float weight[][N_ATTR], int j);
int NO_OF_TRAINING_EXAMPLES=N_SAMPLES;
int NO_OF_INPUTS=N_ATTR-1;
FILE* file_pointer;
FILE* fptr1;
FILE* fptr2;
FILE* output;
//float Vk[sum]={0,0,0,0,0,0,0,0,};
float error(int layer_no, int neuron_no, int  Layer[], float weight[][N_ATTR], float err, float Vj[][1]);
float delta(int layer_no, int neuron_no, int  Layer[], float weight[][N_ATTR], float err, float Vj[][1]);
float Derivative_Phi(float z);

float Forward_Feed(int x, int Layer[], float weight[][N_ATTR], float input[], float Vk[][1], float Y[], float input_storage[][N_ATTR])
{
    int k,i,size;
    float Weighted_Sum;
  float Calculate_Weighted_Sum(float input[],int size,float weight[][N_ATTR], int j);
  float Activation_Function(float Weighted_Sum, int Layer_no, int Layer[]);
  float Check_Output(FILE * file_pointer, float actual_output);
  float error=0;
  float c;
    //int y[Layer[x]+1];
    int j;
    //y[0]=1;
  
    size=Layer[x-1]+1;
    if(x==1)
    {
        input_storage[0][0]=1;
        //input[0]=1;
        for(k=1;k<NO_OF_INPUTS+1;k++)
        {
            fscanf(file_pointer,"%f ",&c);
            //printf("Input %f \n",c);
            input[k]=c;
        	input_storage[0][k]=c;
        }
         k=0;
        j=0;
        Y[0]=1;
        //printf("Weighted Sum is: \n");
       for(i=0;i<Layer[x];i++)
      {
          Weighted_Sum=Calculate_Weighted_Sum(input,size,weight,j);
          //printf("For layer %d and neuron %d: %f \n",x,j,Weighted_Sum);
          Vk[i][0]=Weighted_Sum;
          //printf("Vk[%d]=%f\n",i,Vk[i][0]);
          Y[i+1]=Activation_Function(Weighted_Sum,x,Layer);
          //printf("Output: %f",Y[i+1]);
          j++;
      }
    }
    else
    {
        j=0;
        i=0;
        for(i=1;i<x;i++)
        {
          j=j+Layer[i];
        }
        Y[0]=1;
        i=0;
       //printf("Weighted Sum is: \n");
       for(i=0;i<Layer[x];i++)
        {
            Weighted_Sum=Calculate_Weighted_Sum(input,size,weight,j);
         //   printf("For layer %d and neuron %d: %f \n",x,j,Weighted_Sum);
            Vk[i+j][0]=Weighted_Sum;
           // printf("Vk[%d]=%f\n",i+j,Vk[i+j][0]);
            Y[i+1]=Activation_Function(Weighted_Sum,x,Layer);
            //printf("Output: %f",Y[i+1]);
            
            j++;
        }
      if(Layer[x]==1)
      {
        
        error=Check_Output(file_pointer,Y[1]);
      }
      
    }
   

    for(i=0;i<size;i++)
    {
      if(i<(Layer[x]+1))
      {
        input[i]=Y[i];
      }
      else
      {
        input[i]=0;
      }
    }
    return(error);
    
    
}

float Calculate_Weighted_Sum(float input[], int size, float weight[][N_ATTR], int j)
{
    int i=0;
    float sum=0;
    for( i=0;i<size;i++)
    {
        sum=sum+input[i]*weight[j][i];
    }
    return(sum);
}

float Activation_Function(float Weighted_Sum, int Layer_no, int Layer[])
 {
   float value,x;
   float lambda=LAMBDA;
   if(Layer[Layer_no]!=1)
   {
   		x=(float)(1+exp(-lambda*Weighted_Sum));
   		/*if(x==0)
   		{
   			//printf("Divided By ZERO!");
		   }*/
     //value=Layer[Layer_no]*((float)1/(float)(1+exp(-lambda*Weighted_Sum)));
     value=((float)1/(float)(1+exp(-Weighted_Sum)));
   }
   else
   {
     //printf("Here in last layer!");
     value=floor(1.5*((float)1/(float)(1+exp(-Weighted_Sum))));
   }
   
   return(value);
   
 }

float Check_Output(FILE * file_pointer, float actual_output)
{
  float expected_output=0;
  float error = 0;
  fscanf(file_pointer, "%f ", &expected_output);
  printf("\n");
  fprintf(fptr1,"Expected Output: %f\n", expected_output);
  fprintf(fptr1,"Actual Output: %f\n", actual_output);
  error = expected_output - actual_output;
  return error;
}

void EBP(float weight[][N_ATTR],float Vk[][1], int Layer[], int n_layer, float err,float input_storage[][N_ATTR], int sum)
{
	int i,j,t=0,k,m;
	float correction_in_layer[n_layer],crr;
	correction_in_layer[0]=0;
	float neta=LEARNING_RATE;
	float correction(int layer_no, int n_layer, float err, float Vk[][1], float weight[][N_ATTR],  int y, int x, int Layer[]);
/*	for(m=1;m<n_layer;m++)
	{
		correction_in_layer[m]=correction(m,n_layer,err,Vk,weight,)
	}*/
	t=0;
	i=0;
	for(i=0;i<sum;)
	{
		//crr=correction(t+1,n_layer,err,Vk,weight,i,Layer[t],Layer);
		//float error=layer_error(weight,i,Layer[t],delta(layer_no+1))
		
		if(i<Layer[t+1])
		{
			j=0;
			crr=error(t+1,i,Layer,weight,err,Vk);
			for(j=0;j<Layer[t]+1;j++)
			{
				weight[i][j]=weight[i][j]+(neta*crr*input_storage[t][j]);
			}
			i++;
		}
		else
		{
			t++;
			//crr=error(t+1,i,Layer,weight,err,Vk);
			k=0;
			for(k=0;k<Layer[t+1];k++)
			{
				j=0;
				crr=error(t+1,i+k,Layer,weight,err,Vk);
				//printf("i+k value: %d \n",i+k);
				for(j=0;j<Layer[t]+1;j++)
				{
					//UPDATE WEIGHT!
					weight[i+k][j]=weight[i+k][j]+(neta*crr*input_storage[t][j]);
				}
				
			}
			i=i+k;
		}
	}
	
}

float error(int layer_no, int neuron_no, int  Layer[], float weight[][N_ATTR], float err, float Vj[][1])
{
	float delta_value=0,sum=0,ret_val=0;
	int outgoing_input_count=Layer[layer_no+1],i,j;
	int start,stop;
	if(neuron_no==0||neuron_no==1||neuron_no==2||neuron_no==3)
	{
		start=4;
		stop=start+outgoing_input_count;
	}
	else if(neuron_no==4||neuron_no==5)
	{
		start=6;
		stop=start+outgoing_input_count;
	}
	
	if(Layer[layer_no]==1)
	{
		ret_val=(err*Derivative_Phi(Vj[neuron_no][0]));
	}
	else
	{
		delta_value=delta(layer_no+1,neuron_no,Layer,weight,err,Vj);
		for(i=start;i<stop;i++)
		{
			if(layer_no==1)
			{
				sum=sum+weight[i][neuron_no+1];
			}
			else if(layer_no==2&&neuron_no==4)
			{
				sum=sum+weight[i][0];
			}
			else if(layer_no==2&&neuron_no==5)
			{
				sum=sum+weight[i][1];
			}
			
		}
	
		ret_val=Derivative_Phi(neuron_no)*sum;
		
		
	}
	//printf("\nret_val neuron %d:  %f\n",neuron_no,ret_val);
	return(ret_val);
}
float delta(int layer_no, int neuron_no, int  Layer[], float weight[][N_ATTR], float err, float Vj[][1])
{
	float error_value=0,ret_val=0;
	error_value=error(layer_no,neuron_no,Layer,weight,err,Vj);
	//printf("error valuein delta value %d: %f\n", neuron_no, error_value);
	ret_val=Derivative_Phi(Vj[neuron_no][0])*error_value;
	//printf("delta value %d: %f\n", neuron_no, ret_val);
	return(ret_val);
	
}
/*void Error_Back_Propagation(float weight[][N_ATTR],float Vk[][1], int Layer[], int n_layer, float err,float input_storage[][N_ATTR])
{
  int x,y,i,l,sum=0,j,k;
  float neta=LEARNING_RATE;
  float total_update;
  //float correction(int layer_no, int n_layer, float err, float Vk[][1], float weight[][N_ATTR],  int y, int x);
  
  for(i=n_layer-1;i>0;i--)
  {
    sum=0;
    for(l=1;l<=i;l++)
    {
      sum=sum+Layer[l];
    }
    y=sum-1;
    x=Layer[i-1]+1;
    int m=x;
    total_update=neta*correction(i,n_layer,err,Vk,weight,y,x);
    for(j=y;j>=0;j--)
    {
      //m=x;
      for(k=x-1;k>=0;k--)
      {
        weight[j][k]=weight[j][k]+total_update*input_storage[i-1][k];
        m--;
      }
    }
  }
}*/

float correction(int layer_no, int n_layer, float err, float Vk[][1],float weight[][N_ATTR],  int y, int x, int Layer[])
{
  float ret_val;
  float cons;
  
  float sum=0;
  int i,j;
  if(layer_no==n_layer-1)
  {
    ret_val=Derivative_Phi(Vk[y][0])*err;
  }
  else
  {
    cons=correction(layer_no+1,n_layer,err,Vk,weight,y+Layer[layer_no],Layer[layer_no+1], Layer);
    for(i=y;i<y+x;i++)
    {
      for(j=0;j<Layer[layer_no];j++)
      {
      	sum=sum+cons*weight[i][j];
	  }
      
    }
    ret_val=Derivative_Phi(Vk[y][0])*sum;
  }
  return(ret_val);
}

/*float dell(int layer_no, int n_layer, float err, float Vk[][1],float weight[][N_ATTR],  int y, int x, int Layer[],int neuron_no)
{
	//int temp=y;
  float ret_val;
  float cons;
  float Derivative_Phi(float z);
  float sum=0;
  float dell3=Derivative_Phi(Vk[y][0])*err;
  int i,j;
  if(layer_no==n_layer-1)
  {
    ret_val=Derivative_Phi(Vk[y][0])*err;
  }
  else
  {
    /*cons=dell(layer_no+1,n_layer,err,Vk,weight,y+Layer[layer_no],Layer[layer_no+1], Layer,  neuron_no);
    for(i=y;i<y+x;i++)
    {
      for(j=0;j<Layer[layer_no];j++)
      {
      	sum=sum+cons*weight[i][j];
	  }
      
    }
    if(layer_no==1)
    {
    	dell1=dell2*dell3*(weight[4][y]+[5][y]);
	}
    
    ret_val=Derivative_Phi(Vk[y][0])*sum;
  }
  return(ret_val);
}*/
float Derivative_Phi(float z)
{
	float lambda=LAMBDA;
	return(z*(1-z));
}

int main()
{
    float Forward_Feed(int x, int Layer[], float weight[][N_ATTR], float input[], float Vk[][1], float Y[], float input_storage[][N_ATTR]);
    void Error_Back_Propagation(float weight[][N_ATTR],float Vk[][1], int Layer[], int n_layer, float err, float input_storage[][N_ATTR]);
    void EBP(float weight[][N_ATTR],float Vk[][1], int Layer[], int n_layer, float err,float input_storage[][N_ATTR], int sum);
    
    int i,sum,n_layers,j,k,Epoch=0,iterations,l,m,n;
    n_layers=ceil((float)log(NO_OF_INPUTS+1)/(float)log(2));
    fprintf(fptr1,"Total no of layers : %d\n ",n_layers);
    int Layer[n_layers];
    float input[NO_OF_INPUTS+1], y[NO_OF_INPUTS+1],err,E=0,Eavg=0;
    
  
  
    file_pointer=fopen("Training.txt", "r");
    fptr1=fopen("Training_Output.txt", "w");
    fptr2=fopen("Testing_Output.txt", "w");
    output=fopen("E_avg.txt", "w");
     if (file_pointer == NULL)
    {
        fprintf(stderr, "Cannot open training data file.\n");
        fprintf(stderr, "Check permissions for data file.\n");
        exit(1);
    }
     if (fptr1 == NULL)
    {
        fprintf(stderr, "Cannot open training data file.\n");
        fprintf(stderr, "Check permissions for data file.\n");
        exit(1);
    }
     if (fptr2 == NULL)
    {
        fprintf(stderr, "Cannot open training data file.\n");
        fprintf(stderr, "Check permissions for data file.\n");
        exit(1);
    }
    
    
    i=0;
    for(i=0;i<NO_OF_INPUTS+1;i++)
    {
        if(i==0)
        {
            input[i]=1;
        }
        else
        {
            input[i]=0;
        }
            
    }
    
    
    i=1;
    Layer[0]=NO_OF_INPUTS+1;
    fprintf(fptr1,"Layer %d: %d\n",0,Layer[0]);
    sum=0;
    for(i=1;i<n_layers;i++)
    {
        Layer[i]=ceil((float)log(Layer[i-1])/(float)log(2));
        fprintf(fptr1,"Layer %d: %d\n",i,Layer[i]);
        if(Layer[i]==1)
        {
            n_layers=i+1;
        }
        sum=sum+Layer[i];
    }
    
    
    fprintf(fptr1,"Layers @ oUTPUT %d\n",Layer[n_layers-1]);
  
    if(Layer[n_layers-1]==1)
    {
        fprintf(fptr1,"Layer made perfectly!\n");
    }
  
    fprintf(fptr1,"Total no of neurons: %d",sum);
    
    float  weight[sum][NO_OF_INPUTS+1];
    float Vk[sum][1];
    float input_storage[n_layers][NO_OF_INPUTS+1];
    float Y[NO_OF_INPUTS+1];
    
    i=0;
    j=0;
    for(i=0;i<sum;i++)
    {
        j=0;
        for(j=0;j<NO_OF_INPUTS+1;j++)
        {
            weight[i][j]=((float)(rand()%10)/(float)10);
            if(j==0)
            {
              Vk[i][j]=0;
            }
            if(i<n_layers)
            {
            	input_storage[i][j]=0;
			}
        }
    }
  	/*0.069 0.75 -0.79 0 0.01 0 0.32 -0.34 0 2
0.347 -0.25 0.64 0 -0.51 0.29 36 1.97 0.62 4
0.236 -0.25 0.01 0 0.48 -0.67 29 -1.28 -0.37 1
-0.65 -0.25 -0.35 0 0.01 0.04 40 -0.34 -0.25 1*/
    /*i=0;
    j=0;
    k=0;
    for(k=0;k<n_layers;k++)
    {
      int x=Layer[k];
      for(i=0;i<sum;i++)
      {
        
        for(j=0;j<x;j++)
        {
          weight[i][j]=((float)(rand()%10)/(float)10);
        }
      }
      
    }*/
    
    
    Epoch=0;
    iterations=NO_ITERATIONS;
    float c;
    while(iterations!=0)
    {
        Epoch=Epoch+1;
        E=0;
        fprintf(fptr1,"EPOCH: %d",Epoch);
        i=1;
        k=0;
        for(k=0;k<NO_OF_TRAINING_EXAMPLES;k++)
        {
            for(j=0;j<sum;j++)
            {
              Vk[j][0]=0;
             
            }
            fprintf(fptr1,"TRAINING EXAMPLE %d: \n\n",k);
            for(i=1;i<n_layers;i++)
            {
                
                fprintf(fptr1,"Feed Forwarding %d \n",i);
                err=Forward_Feed(i,Layer,weight,input,Vk,Y,input_storage);
                l=0;
                for(l=0;l<Layer[i]+1;l++)
                {
                	input_storage[i][l]=Y[l];
                	//printf("New input %d: %f",i,input_storage[i][l]);
				}
            }
            if(err==0)
            {
                fprintf(fptr1,"CORRECTLY CLASSIFIED!\n");
            }
            else
            {
              fprintf(fptr1,"MISCLASSIFIED!\n");
              E=E+(0.5*err*err);
              
              EBP(weight, Vk, Layer, n_layers,err, input_storage,sum);
              fprintf(fptr1,"Weights Updated! \n");
               m=0;
	            n=0;
	            for(m=0;m<sum;m++)
	            {
	            	for(n=0;n<NO_OF_INPUTS+1;n++)
	            	{
	            		fprintf(fptr1," %f ",weight[m][n]);
					}
					fprintf(fptr1,"\n");
				}
            }
           
         }
         Eavg=(float)E/(float)NO_OF_TRAINING_EXAMPLES;
         fprintf(output,"%f\n",Eavg);
      iterations--;
      rewind(file_pointer);
    }
    fclose(file_pointer);
    file_pointer=fopen("Testing.txt","r");
    //Testing
    k=0;
    int count=0;
    int NO_OF_TESTING_EXAMPLES=TEST_SAMPLES;
    for(k=0;k<NO_OF_TESTING_EXAMPLES;k++)
        {
        	j=0;
            for(j=0;j<sum;j++)
            {
              Vk[j][0]=0;
             
            }
            fprintf(fptr2,"TESTING EXAMPLE %d: \n\n",k);
            i=1;
            for(i=1;i<n_layers;i++)
            {
                
                fprintf(fptr2,"Feed Forwarding %d \n",i);
                err=Forward_Feed(i,Layer,weight,input,Vk,Y,input_storage);
                l=0;
                for(l=0;l<Layer[i]+1;l++)
                {
                	input_storage[i][l]=Y[l];
                	//printf("New input %d: %f",i,input_storage[i][l]);
				}
            }
            if(err==0)
            {
                fprintf(fptr2,"CORRECTLY CLASSIFIED!\n");
                count++;
            }
            else
            {
              fprintf(fptr2,"MISCLASSIFIED!\n");
    		}
   
    }
    
    fprintf(fptr2,"%d CORRECTLY CLASSIFIED OUT OF %d TESTING EXAMPLES",count, NO_OF_TESTING_EXAMPLES);
    float accuracy=(float)count/(float)NO_OF_TESTING_EXAMPLES;
    fprintf(fptr2,"Accuracy is %f! ",accuracy);
    
}
