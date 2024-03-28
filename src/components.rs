use bevy::window::PrimaryWindow;
use bevy::prelude::*;

pub const PLAYER_SPEED: f32 = 500.0;
pub const BASE_ROAD: (f32, f32, f32)= (0.0, 70.0, 0.0);
pub const BASE_PLAYER: (f32, f32, f32)= (80.0, BASE_ROAD.1 + 60.0, 0.0);
pub const PLAYER_SIZE: f32 = 70.0;


#[derive(Component)]
pub struct Player {
    position: Landscape
}

pub fn spawn_player(mut commands: Commands,
                    asset_server: Res<AssetServer>,
                    landscape_query: Query<&Landscape>) {

    let landscape = landscape_query.get_single().unwrap_or(&Landscape {
        bin_id: 0,
        sv_state: SvState::Ref
    });

    commands.spawn((
        SpriteBundle{
            transform: Transform::from_xyz(BASE_PLAYER.0, BASE_PLAYER.1, BASE_PLAYER.2),
            texture: asset_server.load("sprites/character_maleAdventurer_run0.png"),
            ..Default::default()
        },
        Player {
            position: landscape.to_owned()
        },
        )
    );
}

#[derive(Component)]
pub struct GameCamera {}

pub fn spawn_camera(mut commands: Commands,
                    window_query: Query<&Window, With<PrimaryWindow>>) {
    
    let window = window_query.get_single().unwrap();

    commands.spawn((
        Camera2dBundle {
            transform: Transform::from_xyz(window.width() / 2.0, window.height() / 2.0, BASE_PLAYER.2),
            ..Default::default()
        },
        GameCamera {},
        )
    );
}

pub fn player_movement(keyboard_input: Res<Input<KeyCode>>,
                       mut player_query: Query<(&mut Transform, &Player), With<Player>>,
                       time: Res<Time>,
                       landscape_query: Query<&Landscape, With<Transform>>) {

    if let Ok(mut transform) = player_query.get_single_mut() {
        let mut direction = Vec3::new(0.0, 0.0, BASE_PLAYER.2);

        let speed_offset: f32 = match landscape_query.iter().next().unwrap().sv_state {
            SvState::Ref => 1.0,
            SvState::Amp => 1.0,
            SvState::Del => 0.05,
            SvState::Inv => -1.0
        };

        if keyboard_input.pressed(KeyCode::Right) {
            direction += Vec3::new(1.0, 0.0, BASE_PLAYER.2)
        }
        if keyboard_input.pressed(KeyCode::Left) {
            direction += Vec3::new(-1.0, 0.0, BASE_PLAYER.2)
        }
        if keyboard_input.pressed(KeyCode::Up) {
            direction += Vec3::new(0.0, 1.0, BASE_PLAYER.2)
        }
        if keyboard_input.pressed(KeyCode::Down) {
            direction += Vec3::new(0.0, -1.0, BASE_PLAYER.2)
        }

        if direction.length() > 0.0 {
            direction = direction.normalize();
        }

        transform.0.translation += direction * PLAYER_SPEED * time.delta_seconds();
    }
}

pub fn camera_movement(keyboard_input: Res<Input<KeyCode>>,
                       mut camera_query: Query<&mut Transform, With<GameCamera>>,
                       time: Res<Time>,
                       landscape_query: Query<&Landscape, With<Transform>>) {

    if let Ok(mut transform) = camera_query.get_single_mut() {
        let mut direction = Vec3::ZERO;


        let speed_offset: f32 = match landscape_query.iter().next() {
            Some(landscape) => match landscape.sv_state {
                SvState::Ref => 1.0,
                SvState::Amp => 1.0,
                SvState::Del => 0.05,
                SvState::Inv => -1.0
            },
            None => 2.0
        };

        if keyboard_input.pressed(KeyCode::Right) {
            direction += Vec3::new(1.0, 0.0, 0.0)
        }
        if keyboard_input.pressed(KeyCode::Left) {
            direction += Vec3::new(-1.0, 0.0, 0.0)
        }
        if keyboard_input.pressed(KeyCode::Up) {
            direction += Vec3::new(0.0, 1.0, 0.0)
        }
        if keyboard_input.pressed(KeyCode::Down) {
            direction += Vec3::new(0.0, -1.0, 0.0)
        }

        if direction.length() > 0.0 {
            direction = direction.normalize();
        }

        transform.translation += direction * PLAYER_SPEED * time.delta_seconds() * speed_offset;
    }
}



pub fn confine_player_movement(mut player_query: Query<&mut Transform, With<Player>>,
                               window_query: Query<&Window, With<PrimaryWindow>>) {

    if let Ok(mut player_transform) = player_query.get_single_mut() {
        
        let window = window_query.get_single().unwrap();
        let x_min = BASE_PLAYER.0;
        let x_max = window.width() - PLAYER_SIZE / 2.0;
        let y_min = BASE_PLAYER.1;
        let y_max = y_min + BASE_ROAD.1;

        let mut translation = player_transform.translation;

        //bound the movement
        if translation.x < x_min {
           translation.x = x_min; 
        }
        if translation.x > x_max {
           translation.x = x_max; 
        }
        if translation.y > y_max {
           translation.y = y_max; 
        }
        if translation.y < y_min {
           translation.y = y_min; 
        }

        player_transform.translation = translation;
    }
    
}


pub fn confine_camera_movement(mut camera_query: Query<&mut Transform, With<GameCamera>>,
                               window_query: Query<&Window, With<PrimaryWindow>>) {

    if let Ok(mut camera_transform) = camera_query.get_single_mut() {
        
        let window = window_query.get_single().unwrap();
        let x_min = BASE_PLAYER.0;
        let x_max = window.width() - PLAYER_SIZE / 2.0;
        let y_min = BASE_PLAYER.1;
        let y_max = y_min + BASE_ROAD.1;

        let mut translation = camera_transform.translation;

        //bound the movement
        if translation.x < x_min {
           translation.x = x_min; 
        }
        if translation.x > x_max {
           translation.x = x_max; 
        }
        if translation.y > y_max {
           translation.y = y_max; 
        }
        if translation.y < y_min {
           translation.y = y_min; 
        }

        camera_transform.translation = translation;
    }
    
}

#[derive(Debug, Clone, Copy)]
enum SvState {
    Ref,
    Amp,
    Del,
    Inv 
}


#[derive(Component, Clone, Copy)]
pub struct Landscape {
    bin_id: u32,
    sv_state: SvState
}

pub fn spawn_landscape(mut commands: Commands,
                       asset_server: Res<AssetServer>) {
    

    let offset = 120.0;
    let rounds: u16 = 10;
    let mut sv_type;
    let mut sprite_dir: &str;
    for round in 0..rounds{

        if round < 5 {
            sv_type = SvState::Ref;
            sprite_dir = "sprites/roadTexture_13.png";
        } else {
            sv_type = SvState::Inv;
            sprite_dir = "sprites/roadTexture_85.png";
        }

        // spawning the landscape
        commands.spawn((
                SpriteBundle{
                    transform: Transform::from_xyz(BASE_ROAD.0 + offset * f32::from(round), BASE_ROAD.1, BASE_ROAD.2 - 1.0),
                    texture: asset_server.load(sprite_dir),
                    ..Default::default()
                },
                Landscape {
                    bin_id: u32::from(round),
                    sv_state: sv_type
                },
                )
        );
    }
}
