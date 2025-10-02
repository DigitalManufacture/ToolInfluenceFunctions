function DiffResult=DIFF(VXY,direction,delta)
[N1,N2]=size(VXY);
DiffResult=VXY;
if direction==2 %%along Y
    for i=1:N1
        for j=1:N2
            if i<N1
                DiffResult(i,j)=(VXY(i+1,j)-VXY(i,j))/delta;
            else
                DiffResult(i,j)=(VXY(i,j)-VXY(i-1,j))/delta;
            end
         end
    end
end

if direction==1
    for i=1:N1
        for j=1:N2
            if j<N2
                DiffResult(i,j)=(VXY(i,j+1)-VXY(i,j))/delta;
            else
                DiffResult(i,j)=(VXY(i,j)-VXY(i,j-1))/delta;
            end
         end
    end
end

end