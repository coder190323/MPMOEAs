function MPHH(Global)
% <algorithm> MPHH <A>
% 多种群遗传算法
% div --- 10 --- The number of divisions in each objective
% imgrate --- 0.1 --- The number of individual in each imgrate
% imgrategen --- 20 --- The number of gen in each imgrate
% Na --- 100 --- The capacity of Archive

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    %% Parameter setting
    [div,imgrate,imgrategen,Na] = Global.ParameterSet(10,0.1,20,100);      %获取用户配置参数    
    MP = 2;                                              %子种群的个数（必须为偶数）
    POPNUM = fix(Global.N/MP);          %每个子种群中个体的个数
    IMNUM = fix(POPNUM*imgrate);     %子种群移民个体数量
    Global.N = MP * POPNUM;              %重新定义种群大小
    
    %% Generate random population
    Population = Global.Initialization(Global.N);
    FrontNo = zeros(1,Global.N);            %非支配等级（指标1）
    CrowdDis = zeros(1,Global.N);         %拥挤距离（指标2）
    Pbest  = Population;
    Offspring  = Population;    
    Target  = Population;    
    
    %% 各个种群初始化
    for i=1:MP
        %种群个体
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        %% 不同的子种群使用不同的策略
        if mod(i,2)~=0
           	[~,FrontNo(bird),CrowdDis(bird)] = EnvironmentalSelection(Population(bird),POPNUM);     
        else
            Pbest(bird)      = Population(bird);
            Archive    = UpdateArchive(Population(bird),[],Na,div);
        end
    end
    
    %% Optimization
    while Global.NotTermination(Target) %显示储备集合中的信息
        for i=1:MP
            %子种群个体
            bird = ((i-1)*POPNUM+1):i*POPNUM;
            %% 不同的子种群使用不同的策略
            if mod(i,2)~=0
                MatingPool = TournamentSelection(2,POPNUM,FrontNo(bird),-CrowdDis(bird));
                Offspring(bird)  = Global.Variation(Population(MatingPool),POPNUM,@EAbinary);
                [Population(bird),FrontNo(bird),CrowdDis(bird)] = EnvironmentalSelection([Population(bird),Offspring(bird)],POPNUM);
            else
                REP        = REPSelection(Archive.objs,POPNUM,div);                             %确定全局引导者
                Population(bird) = Global.Variation([Population(bird),Pbest(bird),Archive(REP)],POPNUM,@BBPSO);               %生成下一代种群
                Pbest(bird)      = UpdatePbest(Pbest(bird),Population(bird));                   %更新个体引导者
                Archive    = UpdateArchive(Population(bird),Archive,Na,div);                   %更新储备集
            end
        end
        
        %% 种群迁移
        if mod(Global.evaluated,imgrategen*Global.N)==0
            %% 选择迁移模式
            %     模式1:子种群间好个体相互交换（环）；
            %     模式2:随机交换（环）；
            %     模式3:好的替换坏的（环）；
            %     模式4:好的替换坏的（单向）。
            Population = immigrant(Population,MP,IMNUM,3);
        end
        
        %% 更新全局储备集
        Target    = UpdateArchive(Population,Target,Na,div);                   
    end
end