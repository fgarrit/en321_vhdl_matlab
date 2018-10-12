library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity registre is
port(
   iClock            : in  std_logic;
   iReset            : in  std_logic;
   iEN               : in  std_logic;
   iData             : in  std_logic_vector(7 downto 0);
   oDataValid        : out std_logic;
   oData             : out std_logic_vector(7 downto 0));
end registre;

architecture Behavioral of registre is
begin
process(iClock, iReset)
begin
   if(iReset = '1')   then
      oData <= (others =>'0');
   elsif(iClock'EVENT and iClock = '1')   then
      if (iEN= '1') then
         oDataV <= iData;      
      end if;   
   end if;
end process;   


end Behavioral;