library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity bch_encoder is
port(
   rst            : in  std_logic;
   clk            : in  std_logic;
   i_data         : in  std_logic_vector(3 downto 0);
   i_dv           : in  std_logic;
   o_data         : out std_logic_vector(6 downto 0);
   o_dv           : out std_logic);
end bch_encoder;

architecture Behavioral of bch_encoder is
begin
process(clk, rst)
begin
   if(rst = '1')   then
      o_dv <= '0';
      o_data <= (others =>'0');
   elsif(clk'EVENT and clk = '1')   then
      --if (iEN= '1') then
         o_data(6 downto 3) <= i_data(3 downto 0) ;
         o_data(2) <= i_data(2) xor i_data(1) xor i_data(0) ;
         o_data(1) <= i_data(3) xor i_data(1) xor i_data(0) ;
         o_data(0) <= i_data(3) xor i_data(2) xor i_data(0) ;    
      --end if;   
   end if;
end process;   

end Behavioral;