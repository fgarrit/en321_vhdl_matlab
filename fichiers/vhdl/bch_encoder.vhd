
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity bch_encoder is
port(
   iClock            : in  std_logic;
   iReset            : in  std_logic;
   iEN               : in  std_logic;
   iDataParallel     : in  std_logic_vector(3 downto 0);
   oDataValid        : out std_logic;
   oDataParallel     : out std_logic_vector(6 downto 0));
end bch_encoder;

architecture Behavioral of bch_encoder is
begin
process(iClock, iReset)
begin
   if(iReset = '1')   then
      oDataParallel <= (others =>'0');
   elsif(iClock'EVENT and iClock = '1')   then
      if (iEN= '1') then
         oDataParallel(6 downto 3) <= iDataParallel(3 downto 0) ;
         oDataParallel(2) <= iDataParallel(2) xor iDataParallel(1) xor iDataParallel(0) ;
         oDataParallel(1) <= iDataParallel(3) xor iDataParallel(1) xor iDataParallel(0) ;
         oDataParallel(0) <= iDataParallel(3) xor iDataParallel(2) xor iDataParallel(0) ;    
      end if;   
   end if;
end process;   

end Behavioral;