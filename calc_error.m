function [mae,mse,rmse,mape,error,errorPercent,mre]=calc_error(x1,x2)
%此函数用于计算预测值和实际（期望）值的各项误差指标
%   参数说明
%----函数的输入值-------
%   x1：真实值
%   x2：预测值

%----函数的返回值-------
%   mae：平均绝对误差（是绝对误差的平均值，反映预测值误差的实际情况.）
%   mse：均方误差（是预测值与实际值偏差的平方和与样本总数的比值）
%   rmse：均方误差根（是预测值与实际值偏差的平方和与样本总数的比值的平方根，也就是mse开根号，
%               用来衡量预测值同实际值之间的偏差）
%   mape：平均绝对百分比误差（是预测值与实际值偏差绝对值与实际值的比值，取平均值的结果，可以消除量纲的影响，用于客观的评价偏差）
%   error：误差
%   errorPercent：相对误差
if nargin==2
    if size(x1,2)==1
        x1=x1';  %将列向量转换为行向量
    end
    
    if size(x2,2)==1
        x2=x2';  %将列向量转换为行向量
    end
    
    num=size(x1,2);%统计样本总数
    error=x2-x1;  %计算误差
    errorPercent=abs(error)./x1; %计算每个样本的绝对百分比误差/相对误差
    
    mae=sum(abs(error))/num; %计算平均绝对误差
    mse=sum(error.*error)/num;  %计算均方误差
    rmse=sqrt(mse);     %计算均方误差根
    mape=mean(errorPercent);  %计算平均绝对百分比误差
    mre=max(abs(error)./x1); % 计算最大相对误差，乘以100%
    %结果输出
    
    disp(['平均绝对误差mae为：              ',num2str(mae)])
    disp(['均方误差mse为：                    ',num2str(mse)])
    disp(['均方误差根rmse为：                ',num2str(rmse)])
    disp(['平均相对/绝对百分比误差mape为：   ',num2str(mape*100),' %'])
    disp(['最大相对误差mre为：   ',num2str(mre*100),' %'])
    disp(['个体相对误差小于平均相对误差的体素个数：  ',num2str(sum(sum(errorPercent<mape)))])
    
else
    disp('函数调用方法有误，请检查输入参数的个数')
end

end

