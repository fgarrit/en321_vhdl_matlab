library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity P2S is
generic (width: integer := 7);
port( 
      clk                  : in std_logic;
      reset                : in std_logic;
      load                 : in std_logic;
      par_data             : in std_logic_vector(width-1 downto 0);
      serial_data          : out std_logic;
      serial_data_valid    : out std_logic);
end P2S;

architecture Behavioral of P2S is
signal cmp : integer := 0;
signal stock : std_logic_vector(6 downto 0) := (others => '0'); 
begin

   process(clk, rst)
   begin
     if (reset = '1')   then
         serial_data <= '0';
         serial_data_valid <= '0';

      elsif rising_edge(clk) then

            serial_data_valid <= '0'; -- 0 par défaut

            if load = '1' then
            cmp <= cmp + 1;
            
            if cmp = 7 then
               cmp <= 0;
               serial_data_valid <= '1'; -- 7ème bit reçu, les registres sont prêts et la valeur de sortie est juste

            end if;

            stock <= stock (2 downto 0) & serial_data; -- Registre à décalage

         end if;
      end if;
   end process;

   o_parrallel_data <= stock; -- On relie les registres à la sortie

end Behavioral;