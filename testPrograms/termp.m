
DRNLParams.a=10;
paramChanges={'DRNL.a=2;'}
for idx=1:length(paramChanges)
x=paramChanges{idx};
st=x(1:strfind(x,'.')-1);
fld=x(strfind(x,'.')+1:strfind(x,'=')-1);
x1=eval(['isstruct(' st ')']);
x2=eval(['isfield(' st ',''' fld ''')']);
if ~x1*x2
    disp([' paramChange is faulty: ' x])
end
end