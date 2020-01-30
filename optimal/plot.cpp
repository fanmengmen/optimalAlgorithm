// #include<opencv2/opencv.hpp>
// using namespace cv;
// int main()
// {
// 	//1.从摄像头读入视频
//         Mat cam = imread("motianlun.jpg");
// 		imshow("相机",cam);//显示当前帧图像
// 		waitKey(0);//延时30秒
	
// 	return 0;
// }
#include "plot.hpp"
#define pi 3.1415926

using namespace cv;
using namespace std;
int plot(Mat &M1,double &r, double &c, double &d, double alpha)
{
    int samplePoint = 10;
    int m = 100;
    int n = 100;
    

	  double speed = 2;
    double sampleTime = 0.1;
    vector<double> data_x;
    vector<double> data_y;
    //没有旋转的8字
    double T = 2*r/speed;
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


    // for(int i = 0; i < data_x.size(); i++)
    // {
    //     cout << data_x[i] << "  " << data_y[i] << endl;
    // }


// 带有旋转的8字
// for(double t = 0; t < T; t=t+0.1)
//     {
//       double x = -sqrt(2)/2*r + sqrt(2)/2*speed*t;
//       double y = -sqrt(2)/2*r + sqrt(2)/2*speed*t;
//       data_x.push_back(cos(alpha)*x + sin(alpha)*y + c);
//       data_y.push_back(-sin(alpha)*x + cos(alpha)*y + d);
//     }

//     T = 3*pi*r/2/speed;
//     for(double t = 0; t < T; t=t+0.1)
//     {
//       double theta1 = 3*pi/4 - speed*t/r;
//       double x = sqrt(2)*r + r*cos(theta1);
//       double y = r*sin(theta1);
//       data_x.push_back(x*cos(alpha) + y*sin(alpha) + c);
//       data_y.push_back(-x*sin(alpha) + y*cos(alpha) + d);
//     }

//     T = 2*r/speed;
//     for(double t = 0; t < T; t=t+0.1)
//     {
//       double theta1 = 3*pi/4 - speed*t/r;
//       double x = sqrt(2)/2*r - sqrt(2)/2*speed*t;
//       double y = -sqrt(2)/2*r + sqrt(2)/2*speed*t;
//       data_x.push_back(x*cos(alpha) + y*sin(alpha) + c);
//       data_y.push_back(-x*sin(alpha) + y*cos(alpha) + d);     
//     }

//     T = 3*pi*r/2/speed;
//     for(double t = 0; t < T; t=t+0.1)
//     {
//       double theta2 = pi/4 + speed*t/r;
//       double x = -sqrt(2)*r + r*cos(theta2);
//       double y = r*sin(theta2);
//       data_x.push_back(x*cos(alpha) + y*sin(alpha) + c);
//       data_y.push_back(-x*sin(alpha) + y*cos(alpha) + d);
//     }
    Point center;
    for(int i = 0; i < data_x.size(); i++)
    {
      center.y = (data_y[i] + d)*samplePoint;
      center.x = (data_x[i] + c)*samplePoint;
      circle(M1,center,2,alpha,2);
      
    }
    // waitKey(0);
    
    
}
