library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity S2P is
generic (width: integer := 4);
port(
  clk                  : in std_logic;
  reset                : in std_logic;
  i_data_valid         : in std_logic;
  serial_data          : in std_logic;
  par_data             : out std_logic_vector(width-1 downto 0);
  o_data_valid         : out std_logic);
end S2P;


architecture Behavioral of S2P is
signal cmp      : integer := 0;
signal stock    : std_logic_vector(3 downto 0) := (others => '0'); 
begin

   process(clk, reset)
   begin
     if (reset = '1')   then
         stock <= (others => '0');
         o_data_valid <= '0';

      elsif rising_edge(clk) then

         o_data_valid <= '0'; -- 0 par défaut

         if i_data_valid = '1' then
            cmp <= cmp + 1;
            if cmp = 3 then
               cmp <= 0;
               o_data_valid <= '1'; -- 4ème bit reçu, les registres sont prêts et la valeur de sortie est juste
            end if;
            stock <= stock (2 downto 0) & serial_data; -- Registre à décalage
         end if;
      end if;
   end process;

   par_data <= stock; -- On relie les registres à la sortie

end Behavioral;