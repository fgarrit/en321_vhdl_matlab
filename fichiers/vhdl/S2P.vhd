library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity S2P is 
generic(
    width : integer := 4 );
port( 
   clk            : in std_logic;
   reset          : in std_logic;
   i_data_valid   : in std_logic;
   serial_data    : in std_logic;
   par_data       : out std_logic_vector(3 downto 0);
   o_data_valid   : out std_logic);

end S2P;

architecture Behavioral of S2P is
begin
process(clk, reset)
begin  -- A FAIRE --
   if(reset = '1')   then
      par_data <= (others =>'0');
      o_data_valid <= '0';
   elsif(clk'EVENT and clk = '1')   then
      if (i_data_valid= '1') then
   
      end if;   
   end if;
end process;   

end Behavioral;