----------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.numeric_std.all;

entity codeur_conv is
	port(
		iClock       : in	std_logic;
		iReset       : in	std_logic;
		iEN	    	 : in	std_logic;
		iData        : in	std_logic;
		oDataX       : out std_logic;
		oDataY       : out std_logic
	 );
end codeur_conv;

architecture Behavioral of codeur_conv is
begin  -- A FAIRE --
   process (iClock, iReset) begin
      if (rising_edge(iClock)) then
         if (iReset = '1') then
            oDataX <= '0';
            oDataY <= '0';
         elsif(iEN = '1') then
          
         end if;
      end if;
   end process;
end Behavioral;

