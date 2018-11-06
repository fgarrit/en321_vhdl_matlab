library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity p2s is
    Generic(width             : integer := 7);
    Port ( clk                : in STD_LOGIC;
           reset              : in STD_LOGIC;
           load               : in STD_LOGIC;
           par_data           : in STD_LOGIC_VECTOR (width-1 downto 0);
           serial_data        : out STD_LOGIC;
           serial_data_valid  : out STD_LOGIC);
end p2s;

architecture Behavioral of p2s is

signal Cnt                 : integer range 0 to 9 := 0 ;
signal cnt_horloge         : integer range 0 to 6 ; ---------------
signal start_cnt_h         : std_logic ;
signal load_data           : std_logic;
signal cnt_load            : integer range 0 to 6 := 0 ;
signal data_temp           : std_logic_vector (6 downto 0) := (others => '0');
signal serial_data_t       : std_logic;
signal serial_data_valid_t : std_logic;

begin

Process (reset, clk)
begin
    if reset = '1' then
        data_temp <= (others => '0');
    elsif clk'event and clk = '1' then
        if load = '1' then 
            data_temp <= par_data;
        elsif cnt_horloge < 6 then
            data_temp <= data_temp(5 downto 0) & data_temp(6);
        else 
            data_temp <= data_temp;
        end if;
    end if;
end process;
    
Process (reset, clk)
begin
    if reset = '1' then
        cnt_horloge <= 0;
    elsif clk'event and clk = '1' then
        if load = '1' then 
            cnt_horloge <= 0;
        elsif start_cnt_h = '1' then
            cnt_horloge <= cnt_horloge + 1 ; 
        end if;
    end if;
end process;
    
Process (reset, clk)
begin
    if reset = '1' then
        start_cnt_h <= '0';
    elsif clk'event and clk = '1' then
        if load = '1' then 
            start_cnt_h <= '1';
        elsif cnt_horloge = 6 then
            start_cnt_h <= '0' ; 
        end if;
    end if;
end process;

Process (reset, clk)
begin
    if reset = '1' then
        serial_data <= '0';
        serial_data_valid <= '0';
    elsif clk'event and clk = '1' then
        if start_cnt_h = '1' then 
            serial_data <= data_temp(6);
            serial_data_valid <= '1';
        else
            serial_data_valid <= '0';
        end if;
    end if;
end process;
end Behavioral;