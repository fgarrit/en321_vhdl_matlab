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
signal cmp : integer range 0 to 9 := 0;
signal cmpclk : integer range 0 to 9 := 0;
signal start_cmpclk : std_logic;
signal load_data : std_logic;
signal cmp_load : integer range 0 to 6 := 0;
signal stock : std_logic_vector(6 downto 0) := (others => '0'); 
begin

   process(clk, reset)
   begin
      if (reset = '1')   then
         start_cmpclk <= '0';
      elsif rising_edge(clk) then
            if load = '1' then 
               start_cmpclk <= '1';
            elsif cmp_load = 7 then
               start_cmpclk <= '0'; 
            end if;
      end if;      
   end process;

   process(clk, reset)
   begin
      if (reset = '1')   then
         cmpclk <= 0;
      elsif rising_edge(clk) then
         if start_cmpclk = '1' then
            if cmpclk > 9 then
               cmpclk <= 0;
            else 
               cmpclk <= cmpclk + 1;
            end if;
         else
            cmpclk <= 0;
         end if;      
      end if;   
   end process;

   process(clk, reset) --pause
   begin
      if cmpclk = 10 then
         load_data <= '1';
      else 
         load_data <= '0';           
      end if ;
   end process;

   process (reset, clk)
   begin
      if reset = '1' then
        cmp_load <= 0;
    elsif rising_edge (clk) then
        if load = '1' then
            cmp_load <= 0;
        end if;
        if load_data = '1' then
            if cmp_load > 6 then
                cmp_load <= 0;
            else
                cmp_load <= cmp_load +1;
            end if;
        end if;
    end if;
   end process;

process(reset, clk)
begin
    if reset = '1' then
        stock <= "0000000";
    elsif rising_edge(clk) then
        if load='1' then
            stock <= par_data;
        elsif load_data = '1' then
            stock <= stock( 5 downto 0) & stock(6);
        end if;
    end if;
end process;  

process(load_data, reset)
begin
    if reset = '1' then
        serial_data <= '0';
    elsif load_data ='1' then
        serial_data <= data_temp(6);
    end if;
end process;

serial_data_valid <= load_data;
end Behavioral;