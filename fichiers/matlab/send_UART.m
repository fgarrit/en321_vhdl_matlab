function s = send(vec, size)

s = serial('COM5');
s.Baudrate=115200;
s.StopBits=1;
s.Parity='none';
s.FlowControl='none';
s.TimeOut = 1;
s.InputBufferSize = 100000;
s.OutputBufferSize = 100000;
fopen(s);

for i = 1:size
     fwrite(s,vec(i))
     %pause(0.01) % wait 1ms to let the data arrive to the FPGA
 end
