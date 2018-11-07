library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.numeric_std.all;

entity codeur_conv is
	port(
		iClock       : in	std_logic;
		iReset       : in	std_logic;
		iEN	    	   : in	std_logic;
		iData        : in	std_logic;
		oDataX       : out std_logic;
		oDataY       : out std_logic
	 );
end codeur_conv;

architecture Behavioral of codeur_conv is
signal registre        : std_logic_vector(1 downto 0);
signal s_dataX    : std_logic;
signal s_dataY    : std_logic;
begin 
process(iClock, iReset)
begin
   if(iReset = '1')  then
     registre <= (others=>'0');
   elsif(iClock'EVENT and iClock = '1')   then
      if(iEN = '1')  then
         registre <= iData & registre(1) ;
      end if;        
   end if;
end process;

   oDataX <= iData  xor registre(0);
   oDataY <= registre(1) xor registre(0); 

end Behavioral;