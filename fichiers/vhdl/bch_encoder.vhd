library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity bch_encoder is
port(
   rst            : in  std_logic;
   clk            : in  std_logic;
   i_data               : in  std_logic;
   i_dv            : in  std_logic_vector(3 downto 0);
   o_data        : out std_logic;
   o_dv            : out std_logic_vector(6 downto 0));
end bch_encoder;

architecture Behavioral of bch_encoder is
begin
process(clk, rst)
begin
   if(rst = '1')   then
      o_data <= (others =>'0');
   elsif(clk'EVENT and clk = '1')   then
      --if (iEN= '1') then
         o_data(6 downto 3) <= i_data(3 downto 0) ;
         o_data(2) <= iData(2) xor iData(1) xor iData(0) ;
         o_data(1) <= iData(3) xor iData(1) xor iData(0) ;
         o_data(0) <= iData(3) xor iData(2) xor iData(0) ;    
      --end if;   
   end if;
end process;   

end Behavioral;