function out=savee(i,ZF)
    save(['data/pass',num2str(i),'.dat'], 'ZF');
    out=1;
end