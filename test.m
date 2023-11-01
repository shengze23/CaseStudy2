
clear all;
close all;
% 创建未知变量
syms SS SI IS II IR ID RS a b c d RR

% 创建矩阵和向量

solutions = {};
% 存储所有的解
for SS = 0:0.2:1
    for IS = 0:0.2:1
        for II = 0:0.2:1
            for IR = 0:0.2:1
                for RS = 0:0.2:1
                    for SI = 0:0.2:1
                        for ID = 0:0.2:1
                            for RR= 0:0.2:1
                                for a = 0:0.2:1
                                    for b = 0:0.2:1
                                        for c =0:0.2:1
                                            for d=0:0.2:1
                                if a+b+c+d==1 && SS+SI==1 && IS+II+IR+ID==1 && RS+RR==1
                                    vector(1) = a;
                                    vector(2) = b;
                                    vector(3) = c;
                                    vector(4) = d;
                                    matrix(1, 1) = SS;
                                    matrix(1, 2) = IS;
                                    matrix(1, 3) = RS;
                                    matrix(1, 4) = 0;
                                    matrix(2, 1) = SI;
                                    matrix(2, 2) = II;
                                    matrix(2, 3) = 0;
                                    matrix(2, 4) = 0;
                                    matrix(3, 1) = 0;
                                    matrix(3, 2) = IR;
                                    matrix(3, 3) = RR;
                                    matrix(3, 4) = 0;
                                    matrix(4, 1) = 0;
                                    matrix(4, 2) = ID;
                                    matrix(4, 3) = 0;
                                    matrix(4, 4) = 1;
                                    solutions{end+1} = [SS SI IS II IR ID RS a b c d RR];
                                    for t=2:16
                                        vector=matrix.*vector;
                                        solutions{end+1} = vector;
                                    end
                                else
                                end
                      
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% 显示所有排列组合和解
%for i = 1:2:15    
 %   disp('Matrix * Vector =');
  %  disp(solutions{i + 1});
   % fprintf('\n');
%end
