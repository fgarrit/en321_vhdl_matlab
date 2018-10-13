library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity P2S is
generic (
   width : integer := 4) ;
port (
   clk                  : in std_logic;
   reset                : in std_logic;
   load                 : in std_logic;
   par_data             : in std_logic_vector(width-1 downto 0);
   serial_data          : out std_logic;
   serial_data_valid    : out std_logic);

end P2S;

architecture Behavioral of P2S is
begin
process(clk, reset)
begin  -- A FAIRE --
   if(reset = '1')   then
      serial_data <= '0';
      serial_data_valid <= '0';
   elsif(clk'EVENT and clk = '1')   then
      if (load= '1') then
   
      end if;   
   end if;
end process;   

end Behavioral;