function MPHH(Global)
% <algorithm> MPHH <A>
% ����Ⱥ�Ŵ��㷨
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
    [div,imgrate,imgrategen,Na] = Global.ParameterSet(10,0.1,20,100);      %��ȡ�û����ò���    
    MP = 2;                                              %����Ⱥ�ĸ���������Ϊż����
    POPNUM = fix(Global.N/MP);          %ÿ������Ⱥ�и���ĸ���
    IMNUM = fix(POPNUM*imgrate);     %����Ⱥ�����������
    Global.N = MP * POPNUM;              %���¶�����Ⱥ��С
    
    %% Generate random population
    Population = Global.Initialization(Global.N);
    FrontNo = zeros(1,Global.N);            %��֧��ȼ���ָ��1��
    CrowdDis = zeros(1,Global.N);         %ӵ�����루ָ��2��
    Pbest  = Population;
    Offspring  = Population;    
    Target  = Population;    
    
    %% ������Ⱥ��ʼ��
    for i=1:MP
        %��Ⱥ����
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        %% ��ͬ������Ⱥʹ�ò�ͬ�Ĳ���
        if mod(i,2)~=0
           	[~,FrontNo(bird),CrowdDis(bird)] = EnvironmentalSelection(Population(bird),POPNUM);     
        else
            Pbest(bird)      = Population(bird);
            Archive    = UpdateArchive(Population(bird),[],Na,div);
        end
    end
    
    %% Optimization
    while Global.NotTermination(Target) %��ʾ���������е���Ϣ
        for i=1:MP
            %����Ⱥ����
            bird = ((i-1)*POPNUM+1):i*POPNUM;
            %% ��ͬ������Ⱥʹ�ò�ͬ�Ĳ���
            if mod(i,2)~=0
                MatingPool = TournamentSelection(2,POPNUM,FrontNo(bird),-CrowdDis(bird));
                Offspring(bird)  = Global.Variation(Population(MatingPool),POPNUM,@EAbinary);
                [Population(bird),FrontNo(bird),CrowdDis(bird)] = EnvironmentalSelection([Population(bird),Offspring(bird)],POPNUM);
            else
                REP        = REPSelection(Archive.objs,POPNUM,div);                             %ȷ��ȫ��������
                Population(bird) = Global.Variation([Population(bird),Pbest(bird),Archive(REP)],POPNUM,@BBPSO);               %������һ����Ⱥ
                Pbest(bird)      = UpdatePbest(Pbest(bird),Population(bird));                   %���¸���������
                Archive    = UpdateArchive(Population(bird),Archive,Na,div);                   %���´�����
            end
        end
        
        %% ��ȺǨ��
        if mod(Global.evaluated,imgrategen*Global.N)==0
            %% ѡ��Ǩ��ģʽ
            %     ģʽ1:����Ⱥ��ø����໥������������
            %     ģʽ2:���������������
            %     ģʽ3:�õ��滻���ģ�������
            %     ģʽ4:�õ��滻���ģ����򣩡�
            Population = immigrant(Population,MP,IMNUM,3);
        end
        
        %% ����ȫ�ִ�����
        Target    = UpdateArchive(Population,Target,Na,div);                   
    end
end