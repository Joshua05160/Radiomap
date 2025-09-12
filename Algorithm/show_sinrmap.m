function sinr = show_sinrmap(groundposition, height,pathlossData,N_0)
    [M,N]=size(groundposition);
       xmin_g = -100;
       ymin_g = -160;
       zMin = 10;
       grid_resolution = 5;
       
            for j = 1:N
                groundX = groundposition(1,j);
                groundY = groundposition(2,j);
                groundXIndex = floor((groundX - xmin_g) / grid_resolution) + 1;
                groundYIndex = floor((groundY - ymin_g) / grid_resolution) + 1;
                groundIndex(j) = (groundYIndex-1)*61+groundXIndex;
            end

       airIndex = floor((height - zMin) / grid_resolution)+1;
     for s = 1:N
        for j =1:64
            for i =1:60
               groundIndex_rid_s=groundIndex;
               groundIndex_rid_s(s)=[];
               sinr(i,j,s)=log2(1+10^(-pathlossData(i,j,airIndex,groundIndex(s))/10)/(sum(10.^(-pathlossData(i,j,airIndex,groundIndex_rid_s)/10))+N_0));
           end
       end
     end
end
