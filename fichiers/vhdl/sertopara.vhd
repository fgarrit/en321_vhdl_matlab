library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity sertopara is
port(
   rst               : in  std_logic;
   clk               : in  std_logic;
   i_serial_data     : in  std_logic;
   i_datavalid       : in  std_logic;
   o_parallel_data   : out std_logic_vector(3 downto 0);
   o_dv              : out std_logic);
end sertopara;

architecture Behavioral of sertopara is
signal cmp : integer := 0;
signal stock : std_logic_vector(3 downto 0) := (others => '0'); 
begin

   process(clk, rst)
   begin
     if (iReset = '1')   then
         o_parallel_data <= (others => '0');
         o_dv <= '0';

      elsif rising_edge(clk) then

            o_dv <= '0'; -- 0 par défaut

            if i_datavalid = '1' then
            cmp <= cmp + 1;
            if cmp = 3 then
               cmp <= 0;
               o_dv <= '1'; -- 4ème bit reçu, les registres sont prêts et la valeur de sortie est juste

            end if;

            stock <= stock (2 downto 0) & i_serial_data; -- Registre à décalage

         end if;
      end if;
   end process;

   o_parrallel_data <= stock; -- On relie les registres à la sortie

end Behavioral;