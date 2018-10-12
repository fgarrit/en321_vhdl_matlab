function vec = recv_UART(s, size)

  vec = zeros(1,size);
  
  for i = 1:size
       vec(i) = fread(s,1);
       %pause(0.01) % wait 1ms to let the data arrive to the FPGA
  end

fclose(s);
delete(s)
clear s