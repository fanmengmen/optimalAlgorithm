#include<iostream>
#include<math.h>
#include<vector>
#define pi 3.1415926

using namespace std;

int main()
{
    double speed = 2;
    double r = 15;
    double sampleTime = 0.1;
    vector<double> data_x;
    vector<double> data_y;
    
    double T = 2*r/speed;

    //No.1 直线加圆弧轨迹构建
    for(double t = 0; t < T; t=t+0.1)
    {
      data_x.push_back(-sqrt(2)/2*r + sqrt(2)/2*speed*t);
      data_y.push_back(-sqrt(2)/2*r + sqrt(2)/2*speed*t);
    }

    T = 3*pi*r/2/speed;
    for(double t = 0; t < T; t=t+0.1)
    {
      double theta1 = 3*pi/4 - speed*t/r;
      data_x.push_back(sqrt(2)*r + r*cos(theta1));
      data_y.push_back(r*sin(theta1));
    }

    T = 2*r/speed;
     for(double t = 0; t < T; t=t+0.1)
    {
      double theta1 = 3*pi/4 - speed*t/r;
      data_x.push_back(sqrt(2)/2*r - sqrt(2)/2*speed*t);
      data_y.push_back(-sqrt(2)/2*r + sqrt(2)/2*speed*t);
      
    }

    T = 3*pi*r/2/speed;
    for(double t = 0; t < T; t=t+0.1)
    {
      double theta2 = pi/4 + speed*t/r;
      data_x.push_back(-sqrt(2)*r + r*cos(theta2));
      data_y.push_back(r*sin(theta2));      
    }

    for(int i = 0; i < data_x.size(); i++)
    {
        cout << data_x[i] << "  " << data_y[i] << endl;

    }



  //No.2 双扭线形‘8’字
  // int kNumObservations = 100; //100对点   
  // double data[2*kNumObservations];
  // double theta;
  // double ro;
  // for(int i = 0; i < kNumObservations; i++)
    // {
    //     theta = -pi/4 + 0.01*i + 0.01;
    //     ro = sqrt(cos(2*theta));
    //     data[2*i] = a*ro*cos(theta) + c;
    //     data[2*i+1] = b*ro*sin(theta) + d;
    // }
    


}